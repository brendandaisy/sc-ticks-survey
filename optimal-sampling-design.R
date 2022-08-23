library(INLA)
library(tidyverse)
library(sf)
library(MASS, exclude=c("select"))
library(furrr)
library(gridExtra)

source("other-helpers.R")

parks_model_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data()

prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.5)))

pres_formula <- pres ~ 0 + tick_class +  
     tick_class:land_cover + tick_class:elevation +
    tick_class:tree_canopy + tick_class:jan_min_temp + tick_class:max_temp + 
    tick_class:precipitation + tick_class:mean_rh +
    f(site, model="iid", hyper=prec_pri) +
    # f(month, model="ar1", hyper=prec_pri) + # was eliminated when normals were used
    f(id, model="iid3d", n=nrow(pres), hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3)))))

post_fit <- fit_model(pres_formula, parks_model_data, response="pres", sampling=TRUE)

# Expected risk and variance along grid-------------------------------------------

covar_grid <- st_read("geo-files/covar-grid-16km.shp") |> 
    rename(
        land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
        precipitation=prcpttn, jan_min_temp=jn_mn_t
    ) |> 
    mutate(land_cover=land_cover_labels(land_cover)) |> 
    select(-min_temp)

visited_lc <- levels(parks_model_data$land_cover)
    
covar_grid <- covar_grid |> 
    # group_by(month) |> 
    # slice_sample(n=100) |> 
    # ungroup() |> 
    mutate(
        site=str_c("g", 1:n()),
        data=list(tibble(tick_class=unique(pres$tick_class), pres=NA, tcnum=1:3))
    ) |> 
    unnest(c(data))

pred_loc <- filter(covar_grid, as.character(land_cover) %in% visited_lc)
na_loc <- filter(covar_grid, !(as.character(land_cover) %in% visited_lc))

pred_model_data <- parks_model_data |>
    bind_rows(pred_loc) |> 
    arrange(tcnum) |> 
    mutate(
        across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)), #TODO: decide how scaling is gonna be and remember its diff here!!!
        id=1:n()
    )

pres_formula_grid <- pres ~ 0 + tick_class +  
    tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation +
    tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh +
    # f(site, model="iid", hyper=prec_pri) +
    # f(month, model="ar1", hyper=prec_pri) +
    f(id, model="iid3d", n=nrow(pred_model_data), hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3)))))

fit_grid <- fit_model(pres_formula_grid, pred_model_data, response="pres", sampling=FALSE)

## Get posterior predictive summaries for the prediction grid
post_pred_grid <- pred_model_data |> 
    filter(is.na(pres)) |> 
    mutate(
        mean_pres=fit_grid$summary.fitted.values$mean[which(is.na(pred_model_data$pres))],
        sd_pres=fit_grid$summary.fitted.values$sd[which(is.na(pred_model_data$pres))],
        var_eta=fit_grid$summary.linear.predictor$sd[which(is.na(pred_model_data$pres))]^2
    ) |> 
    bind_rows(na_loc) # add NA locations back in

## Prediction maps based on the observed data
pred_map_theme <- function() {
    theme_bw() +
    theme(
        axis.text=element_blank(), axis.ticks=element_blank(),
        panel.spacing=unit(0, "mm"),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
        legend.position="bottom",
        legend.box.margin=unit(c(0, 0, 0, 0), "mm")
    )
}

pred_map_dat <- post_pred_grid |> 
    mutate(
        month=lubridate::month(month, label=TRUE), # month to english
        tick_class=str_replace(tick_class, " ", "\n")
    )

p1 <- ggplot(pred_map_dat, aes(col=mean_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c() +
    # new_scale_color() +
    # geom_sf(aes(col=as.factor(pres)), data=pres, size=1.5) +
    # scale_color_manual(values=c(`0`="pink", `1`="red")) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Prob. tick present") +
    pred_map_theme()

p2 <- ggplot(pred_map_dat, aes(col=sd_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c(option="magma") +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Prediction std. dev.") +
    pred_map_theme()

ggsave("figs/presence-prediction-map.pdf", grid.arrange(p1, p2, nrow=2), width=12, height=8)

##

spatial_avg <- function(post, var=var_eta) {
    ret <- map_dbl(c(".*", "Amb", "Derm", "Ixo"), ~{
        post |> 
            filter(str_detect(tick_class, .x)) |> 
            pull({{var}}) |> 
            mean(na.rm=TRUE) # mean since the number of prediction locs not constant
    })
    
    names(ret) <- c("All", "Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp.")
    return(ret)
}

spatial_avg(post_pred_grid, var_eta)

##

mdat <- pred_model_data
#TODO: not necessary since this isn't adaptive (see comment below)
mdat$vis <- ifelse(is.na(pred_model_data$pres), TRUE, FALSE)
unvis <- filter(mdat, vis)$site |> 
    unique() |> 
    map(~which(mdat$site == .x))

extract_util <- function(score_list) {
    fext <- function(tick_class) map_dbl(score_list, tick_class)
    map_dbl(c("All", "Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp."), ~mean(fext(.x)))
}

#TODO: not considering dealing with choosing 2nd+ locations, but needing to still "avg over" chosen loc so far (don't see that data)
loc_utility <- function(obs_df, obs_fit, unvis_loc, loc_id, n=10) {
    inew <- unvis_loc[[loc_id]]
    score_tc <- future_map(1:n, ~{
        obs_df$pres[inew] <- rbinom(3, rep(1, 3), obs_fit$summary.fitted.values$mean[inew])
        ft <- fit_model(pres_formula_grid, obs_df, response="pres", sampling=FALSE)
        obs_df |> 
            mutate(var_eta=ft$summary.linear.predictor$sd^2) |> 
            filter(is.na(pres)) |> 
            spatial_avg()
    }, .options=furrr_options(seed=TRUE))
    u_tc <- extract_util(score_tc)
    return(u_tc)
}

plan(multisession, workers=1)
loc_utility(pred_model_data, fit_grid, unvis, 1, n=8)

post_pred_grid <- mdat |> 
    mutate(
        mean_pres=fit$summary.fitted.values$mean,
        sd_pres=fit$summary.fitted.values$sd,
        var_eta=fit$summary.linear.predictor$sd^2
    ) |> 
    filter(is.na(pres))

pred_map_dat <- post_pred_grid |> 
    mutate(
        month=lubridate::month(month, label=TRUE), # month to english
        tick_class=str_replace(tick_class, " ", "\n")
    )

p1 <- ggplot(pred_map_dat, aes(col=mean_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c() +
    # new_scale_color() +
    # geom_sf(aes(col=as.factor(pres)), data=pres, size=1.5) +
    # scale_color_manual(values=c(`0`="pink", `1`="red")) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Prob. tick present") +
    pred_map_theme()

p2 <- ggplot(pred_map_dat, aes(col=sd_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c(option="magma") +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Prediction std. dev.") +
    pred_map_theme()

grid.arrange(p1, p2, nrow=2)

##

library(brms)

tmp <- brm(
    pres ~ 0 + tick_class + tick_class:land_cover + tick_class:elevation + 
        tick_class:tree_canopy + tick_class:jan_min_temp + tick_class:max_temp + 
        tick_class:precipitation + tick_class:mean_rh,
    data=parks_model_data,
    family=bernoulli(),
    prior=c(
        set_prior("normal(0, 1.73)", class="b")
    )
)

##

post_samples <- inla.posterior.sample(1000, post_fit)

presence_predict <- function(..., new_loc=NA) {
    cov <- matrix(c(
        1 / theta[4], theta[7] / sqrt(theta[4]*theta[5]), theta[8] / sqrt(theta[4]*theta[6]),
        theta[7] / sqrt(theta[4]*theta[5]), 1 / theta[5], theta[9] / sqrt(theta[5]*theta[6]),
        theta[8] / sqrt(theta[4]*theta[6]), theta[9] / sqrt(theta[5]*theta[6]), 1 / theta[6]
    ), nrow=3, byrow=TRUE)

    e <- MASS::mvrnorm(1, rep(0, 3), cov)
    s <- rnorm(1, 0, sqrt(1 / theta[1]))
    fx <- get_fixed_effects(..., new_loc=new_loc)

    # notice month-2 is the INDEX of the correct month effect
    eta <- fx + month[new_loc$month-2] + s + e
    r <- inla.link.logit(eta, inverse=TRUE)
    ret <- rbinom(3, rep(1, 3), r)
    names(ret) <- c("Amblyomma americanum", "Demacentor variabilis", "Ixodes spp.")
    return(ret)
}

presence_predict_check <- function(...) {
    r <- inla.link.logit(Predictor, inverse=TRUE)
    m <- length(r) %/% 3
    ret <- rbinom(r, 1, r)
    names(ret) <- c(
        rep("Amblyomma americanum", m), 
        rep("Demacentor variabilis", m),
        rep("Ixodes spp.", m)
    )
    return(ret)
}

# Test sampling with posterior prediction-----------------------------------------

covar_list <- pres |> 
    st_drop_geometry() |> 
    distinct(date, site, .keep_all=TRUE) |> 
    rowwise() |> 
    group_map(~{
        row <- as.list(.x[1,c(2, 5:10)])
        names(row) <- colnames(.x)[c(2, 5:10)]
        row
    })

plan(multisession, workers=4)

inla.posterior.sample.eval(presence_predict_check, post_samples) |> 
    as_tibble(.name_repair="minimal") |> 
    bind_cols(pres) |> 
    rowwise() |> 
    mutate(acc=sum(c_across(contains("sample")) == pres) / length(post_samples)) |> 
    ungroup() |> 
    select(-contains("sample")) |> 
    ggplot(aes(interaction(date, site), acc, col=site)) +
    geom_point(alpha=0.5) +
    facet_wrap(~tick_class, ncol=1) +
    theme(legend.position="none")


pres_predict_list <- future_map(
    covar_list, 
    ~inla.posterior.sample.eval(presence_predict, post_samples, new_loc=.x),
    .options=furrr_options(seed=NULL)
)

tru_pres <- pres |> 
    distinct(date, site, tick_class, pres) |> 
    arrange(date, site)

tru_pres |> 
    bind_cols(
        map_dfr(
            pres_predict_list, ~{
                # enframe(rowMeans(.x))
                as_tibble(.x, .name_repair="minimal")
        })
    ) |> 
    rowwise() |> 
    mutate(acc=sum(c_across(contains("sample")) == pres) / length(pres_predict_list)) |> 
    ungroup() |> 
    select(-contains("sample")) |> 
    ggplot(aes(interaction(date, site), acc, col=site)) +
    geom_point(alpha=0.5) +
    facet_wrap(~tick_class, ncol=1) +
    theme(legend.position="none")

# ^Hard to draw conclusions from this other than did better than random, and probably did will with D. var. by tending to guess
# negative. Weird that 0.6 is such a ceiling
