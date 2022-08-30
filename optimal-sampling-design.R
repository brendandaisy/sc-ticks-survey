library(INLA)
library(tidyverse)
library(sf)
library(MASS, exclude=c("select"))
library(furrr)
library(gridExtra)
library(DiceOptim)

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
    f(id, model="iid3d", n=nrow(parks_model_data), hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3)))))

post_fit <- fit_model(pres_formula, parks_model_data, response="pres", sampling=TRUE)

# Expected risk and variance along grid-------------------------------------------

covar_grid <- st_read("geo-files/covar-grid.shp") |> 
    rename(
        land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
        precipitation=prcpttn, jan_min_temp=jn_mn_t
    ) |> 
    mutate(land_cover=land_cover_labels(land_cover)) |> 
    select(-min_temp)

visited_lc <- levels(parks_model_data$land_cover)

covar_grid <- covar_grid |> 
    group_by(geo=as.character(geometry)) |> 
    mutate(site=str_c("g", cur_group_id())) |> 
    ungroup() |> 
    select(-geo) |> 
    mutate(
        data=list(tibble(tick_class=unique(parks_model_data$tick_class), pres=NA, tcnum=1:3))
    ) |> 
    unnest(c(data))

ggplot(covar_grid, aes(col=precipitation)) +
    geom_sf(size=.2) +
    theme_bw() +
    facet_wrap(~lubridate::month(month, label=TRUE), nrow=1) +
    labs(col="Monthly precipitation") +
    scale_color_viridis() +
    theme(legend.position="bottom", axis.text=element_blank(), axis.ticks=element_blank())

ggsave("figs/precip-monthly.pdf", width=8.5, height=2.2)

pred_loc <- filter(covar_grid, as.character(land_cover) %in% visited_lc)
na_loc <- filter(covar_grid, !(as.character(land_cover) %in% visited_lc))

pred_model_data <- parks_model_data |>
    bind_rows(pred_loc) |> 
    arrange(tcnum) |> 
    mutate(
        across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)), #TODO: decide how scaling is gonna be and remember its diff here!!!
        id=1:n()
    )

get_formula <- function(n) {
    pres ~ 0 + tick_class +  
        tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation +
        tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh +
        f(site, model="iid", hyper=prec_pri) +
        # f(month, model="ar1", hyper=prec_pri) +
        f(id, model="iid3d", n=n, hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3)))))
}

pres_formula_grid <- get_formula(nrow(pred_model_data))
fit_grid <- fit_model(pres_formula_grid, pred_model_data, response="pres", sampling=FALSE)

## Get posterior predictive summaries for the prediction grid
post_pred <- pred_model_data |> 
    # filter(is.na(pres)) |> 
    mutate(
        mean_pres=fit_grid$summary.fitted.values$mean,
        sd_pres=fit_grid$summary.fitted.values$sd,
        var_eta=fit_grid$summary.linear.predictor$sd^2
    )
    # bind_rows(na_loc) # add NA locations back in

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

pred_map_dat <- . |> 
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

###

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

# extract_util <- function(score_list) {
#     fext <- function(tick_class) map_dbl(score_list, tick_class)
#     tc <- c("All", "Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp.")
#     ret <- map_dbl(tc, ~-mean(fext(.x)))
#     names(ret) <- c("All", "Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp.")
#     return(ret)
# }

# sample_score <- function(grid_df, park_df, new_loc_idx, m) {
#     loc_df <- grid_df[new_loc_idx,]
#     loc_df$pres <- rbinom(3, rep(1, 3), loc_df$mean_pres) #TODO: what is relationship between marg risk and the classic nested integral? Confirm same
#     
#     sub_df <- grid_df[-new_loc_idx, ] |> 
#         group_by(month, tick_class) |> 
#         slice_sample(n=m) |> 
#         ungroup() |> 
#         bind_rows(loc_df, park_df) |> 
#         arrange(tcnum) |> 
#         mutate(id=1:n())
#     
#     f <- get_formula(nrow(sub_df))
#     ft <- fit_model(f, sub_df, response="pres", sampling=FALSE)
#     
#     sub_df |> 
#         mutate(var_eta=ft$summary.linear.predictor$sd^2) |> # get the new posterior variance
#         filter(is.na(pres)) |> 
#         spatial_avg()
# }

sample_score <- function(pred_df, park_df, new_loc, m) {
    new_loc$pres <- rbinom(3, rep(1, 3), new_loc$mean_pres) #TODO: what is relationship between marg risk and the classic nested integral? Confirm same
    
    sub_df <- pred_df |> 
        group_by(t, tick_class) |> 
        slice_sample(n=m) |> 
        ungroup() |> 
        bind_rows(new_loc, park_df) |> 
        arrange(tcnum) |> 
        mutate(id=1:n())
    
    f <- get_formula(nrow(sub_df))
    ft <- fit_model(f, sub_df, response="pres", sampling=FALSE)
    
    sub_df |> 
        mutate(var_eta=ft$summary.linear.predictor$sd^2) |> # get the new posterior variance
        filter(is.na(pres)) |> 
        pull(var_eta) |> 
        mean()
}

#TODO: not considering dealing with choosing 2nd+ locations, but needing to still "avg over" chosen loc so far (don't see that data)
# loc_utility <- function(grid_df, park_df, unvis_loc, loc_id, n=10, m=5) {
#     inew <- unvis_loc[[loc_id]]
#     loc_info <- select(grid_df[inew[1],], -tick_class)
#     score_tc <- map(1:n, ~sample_score(grid_df, park_df, inew, m))
#     u_tc <- extract_util(score_tc)
#     return(tibble(enframe(u_tc, "tick_class", "utility"), loc_info))
# }

loc_utility <- function(pred_df, park_df, x, n=10, m=5) {
    d <- fields::rdist(x, as.matrix(select(pred_df, X, Y, t)))
    i_nloc <- which(d == min(d))
    nloc <- pred_df[i_nloc,]
    pred_df <- pred_df[-i_nloc,] # nloc now removed from prediction grid for all reps
    score <- map_dbl(1:n, ~sample_score(pred_df, park_df, nloc, m))
    return(tibble(utility=-mean(score), select(nloc, -tick_class, -tcnum)[1,]))
}

init_utils <- function(pred_df, park_df, n_init=10, n=10, m=5) {
    eval_pts <- pred_df |> 
        slice_sample(n=n_init) |> #DANGER: assumes (correctly ATM) that all points in space have equal weight (each grid point appears 10 t X 3 tick_class)
        select(X, Y, t) |> 
        as.matrix()
        
    map_dfr(
        1:n_init, 
        ~loc_utility(pred_df, park_df, t(eval_pts[.x,]), n, m)
        # .options=furrr_options(seed=TRUE)
    )
}

library(fields)

# post_pred_grid <- filter(post_pred, is.na(pres))
# post_pred_parks <- filter(post_pred, !is.na(pres))

# list of indices for each tick_class, for each grid location
# unvis <- post_pred_grid$site |> 
#     unique() |> 
#     map(~which(post_pred_grid$site == .x))
# 
# tmp <- loc_utility(post_pred_grid, post_pred_parks, unvis, 2022, n=10, m=1)
# 
# plan(multisession, workers=2)
# plan(sequential)

sp_norm <- function(df) {
    df_norm <- df |> 
        st_normalize() |> 
        mutate(month=(month-3) / 9)
}

#Note have to normalize whole prediction space, since using just initial utils as before would mean
# Alg will not check outside bounds of initial utils (which could be much smaller than full space)
pp_norm <- sp_norm(post_pred)

pp_norm <- pp_norm |> 
    st_drop_geometry() |> 
    bind_cols(st_coordinates(pp_norm)) |> 
    rename(t=month)

pp_grid <- filter(pp_norm, is.na(pres))
pp_parks <- filter(pp_norm, !is.na(pres))

loc_utility(pp_grid, pp_parks, res$par, n=1, m=10)

utils0 <- init_utils(pp_grid, pp_parks, n_init=8, n=10, m=20)

krig <- km(
    ~1, select(utils0, X, Y, t), utils0$utility, scaling=TRUE,
    noise.var=rep(0.1, nrow(utils0)), lower=rep(0.1, 3), upper=rep(1, 3)
)

ggplot(utils0, aes(X, Y, col=utility)) +
    geom_point() +
    facet_wrap(~as.factor(t), nrow=2)

res <- max_EQI(krig, new.noise.var=0.1, type="UK", lower=rep(0, 3), upper=rep(1, 3), control=list(max.generations=20))

next_u_eval <- loc_utility(pp_grid, pp_parks, res$par, n=10, m=20)
ucur <- bind_rows(utils0, next_u_eval)

krig2 <- km(
    ~1, select(ucur, X, Y, t), ucur$utility, scaling=TRUE,
    noise.var=rep(0.1, nrow(ucur)), lower=rep(0.1, 3), upper=rep(1, 3)
)
res <- max_EQI(krig2, new.noise.var=0.1, type="UK", lower=rep(0, 3), upper=rep(1, 3), control=list(max.generations=20))

u_eval2 <- loc_utility(pp_grid, pp_parks, res$par, n=10, m=20)
ucur <- bind_rows(ucur, u_eval2)

krig3 <- km(
    ~1, select(ucur, X, Y, t), ucur$utility, scaling=TRUE,
    noise.var=rep(0.1, nrow(ucur)), lower=rep(0.1, 3), upper=rep(1, 3)
)
res <- max_EQI(krig3, new.noise.var=0.3, type="UK", lower=rep(0, 3), upper=rep(1, 3))

dom <- distinct(pp_grid, X, Y, t)
krig_pred <- predict(krig3, dom, type="UK", bias.correct=TRUE)



ss <- length(unique(dom$X))
dom_arr <- array(c(dom$X, dom$Y, dom$t), dim=c(ss, ss, 10))

dom_arr[,,1]

library(akima)
df <- mutate(dom, kmean=krig_pred$mean)
fld <- map(unique(df$t), ~with(filter(df, t == .x), interp(x = X, y = Y, z = kmean)))

map_dfr(fld, ~)

library(ggnewscale)

ggplot(slice_sample(df, n=10000), aes(X, Y, col=kmean)) +
    geom_point(size=1.5, alpha=0.2) +
    facet_wrap(~as.factor(round(t, 2)), nrow=2) +
    scale_color_continuous(guide="none") +
    new_scale_color() +
    geom_point(aes(X, Y, col=utility), data=ucur, size=2, inherit.aes=FALSE) +
    geom_point(aes(X, Y), data=ucur, col="white", size=2.05, shape=1, inherit.aes=FALSE) +
    geom_point(aes(X, Y), data=ucur[9,], col="green", shape=1, size=4.5, inherit.aes=FALSE) +
    geom_point(aes(X, Y), data=ucur[10,], col="green", shape=1, size=4.5, inherit.aes=FALSE) +
    scale_color_viridis_c(option="magma") +
    theme_bw() +
    theme(axis.ticks=element_blank(), axis.text=element_blank())

ggsave("figs/bo-prog.pdf", width=7.2, height=3.6)

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
