library(INLA)
library(tidyverse)
library(sf)
library(MASS, exclude=c("select"))
library(furrr)
library(gridExtra)
library(DiceOptim)
library(fields)
library(hash)

source("other-helpers.R")

parks_model_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

# prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.5)))

# Expected risk and variance along grid-------------------------------------------

covar_grid <- st_read("geo-files/covar-grid.shp") |> 
    rename(
        land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
        precipitation=prcpttn, jan_min_temp=jn_mn_t
    ) |> 
    mutate(land_cover=land_cover_labels(land_cover)) |> 
    select(-min_temp)

# visited_lc <- levels(parks_model_data$land_cover)

covar_grid <- covar_grid |> 
    group_by(geo=as.character(geometry)) |> 
    mutate(site=str_c("g", cur_group_id())) |> 
    ungroup() |> 
    select(-geo) |> 
    mutate(data=list(tibble(tick_class=unique(parks_model_data$tick_class), pres=NA))) |> 
    unnest(c(data))

pred_model_data <- parks_model_data |>
    bind_rows(covar_grid) |>
    arrange(date, site) |>
    mutate(across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)))

f_no_corr <- pres ~ 0 + tick_class +  
    tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation +
    tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh

fit_grid <- fit_model(f_no_corr, st_drop_geometry(pred_model_data), response="pres", fx_prec=0.05)

## Get posterior predictive summaries for the prediction grid
post_pred <- pred_model_data |> 
    mutate(
        mean_pres=fit_grid$summary.fitted.values$mean,
        sd_pres=fit_grid$summary.fitted.values$sd,
        var_eta=fit_grid$summary.linear.predictor$sd^2
    )

# currently not used (only "All" is currently calculated), which means util would have to be recomputed if want
# to do scoring for specific species. This only effects the initial util0 though, since no way to save much work
# after the decisions diverge for partic. species
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

sp_norm <- function(df) {
    df_norm <- df |> 
        st_normalize() |> 
        mutate(month=(month-3) / 9)
}

sample_score <- function(pred_loc, known_loc, samp_loc, m=NULL) {
    # prepare subset of data to make predictions cond. on park_df and new y
    # NOTE: some land_cover classes will not be in this set, but thats OK,
    # because the unshrunk priors from other LCs would not be taken into account
    # when using this subset for sp. avg. variance, by definition
    if (is.null(m)) {
        sub_df <- bind_rows(pred_loc, samp_loc, known_loc)
    }
    else {
        sub_df <- pred_loc |> 
            group_by(t, tick_class) |> 
            slice_sample(n=m) |> 
            ungroup() |> 
            bind_rows(new_loc, known_loc)
    }

    ft <- fit_model(f_no_corr, sub_df, response="pres", fx_prec=0.05)
    
    sub_df |> 
        mutate(var_eta=ft$summary.linear.predictor$sd^2) |> # get the new posterior variance
        filter(is.na(pres)) |> 
        pull(var_eta) |> 
        mean()
}

loc_utility <- function(
        x, pred_loc, known_loc, n=10, m=NULL, 
        by=c("date", "site"), silent=FALSE, list_col=FALSE
    ) {
    # from locations for prediction, find row closest to proposed point x:
    # d <- fields::rdist(x, as.matrix(select(pred_loc, X, Y, t)))
    # i_nloc <- which(d == min(d))
    # add proposed loc to new locs in design, and remove it from prediction set for all reps
    # proposed_loc <- pred_loc[i_nloc,]
    x_loc <- semi_join(pred_loc, x, copy=TRUE, by=by)
    pred_loc <- anti_join(pred_loc, x, copy=TRUE, by=by)
    samp_cache <- hash()
    scores <- map_dbl(1:n, ~{
        # sample new y | new_loc:
        x_loc$pres <- rbinom(nrow(x_loc), rep(1, nrow(x_loc)), x_loc$mean_pres)
        samp_str <- str_c(x_loc$pres, collapse="")
        if (has.key(samp_str, samp_cache)) {
            s <- samp_cache[[samp_str]]
        } else {
            s <- sample_score(pred_loc, known_loc, x_loc, m)
            if (is.null(m))
                samp_cache[[samp_str]] <- s
        }
        if (!silent)
            print(paste0("Score from sample ", .x, "/", n, ": ", s))
        return(s)
    })
    if (list_col)
        return(tibble_row(utility=mean(scores), design=list(x)))
    return(tibble(utility=mean(scores), x))
}

save_util_res <- function(utils, exper=c("random", "bayes_opt"), n, m=NULL, append=TRUE) {
    if (!is.null(m))
        return ("Are you sure? Only should save with full grid")
    file <- paste0("util-exper/", exper, "-n=", n, ".rds")
    if (file.exists(file))
        res <- readRDS(file)
    else
        res <- tibble()
    saveRDS(bind_rows(res, utils), file)
}

# while only necessary for BO (not even RN since scaling anyways), can do this now
pp_norm <- sp_norm(post_pred)

pp_norm <- pp_norm |> 
    st_drop_geometry() |> 
    bind_cols(st_coordinates(pp_norm)) |> 
    rename(t=month)

pp_grid <- filter(pp_norm, is.na(pres))
pp_parks <- filter(pp_norm, !is.na(pres))

# Entire utility surface for |d|=1------------------------------------------------

all_loc <- pp_grid |> 
    distinct(site, date) |> 
    rowwise() |> 
    group_split()

single_d_grid <- map_dfr(all_loc, ~loc_utility(.x, pp_grid, pp_parks, n=20, m=NULL))

# Random selection, baseline design strategy--------------------------------------

B <- 6
reps <- 5
expg <- expand_grid(B=seq(11, 16, 5), rep=1:reps)

u_random <- pmap_dfr(expg, ~{
    d <- slice_sample(pp_grid, n=..1) |> select(date, site)
    mutate(loc_utility(d, pp_grid, pp_parks, n=20, m=NULL, list_col=TRUE), B=..1)
})

# u_random <- u_random |> 
#     mutate(B=c(rep(1, 5), rep(6, 30)), run_on=lubridate::today()) |>
#     select(-rep) |> 
#     nest(design=c(date, site)) |> 
#     select(utility, design, B, run_on)

save_util_res(u_random, "random", n=20)

# Simple variance strategy:-------------------------------------------------------
# Pick the B new locations with highest variance based on data seen---------------

B <- 100

pp_grid_sort <- pp_grid |> 
    group_by(date, site) |> 
    summarise(var_eta=mean(var_eta), .groups="drop") |> 
    arrange(desc(var_eta)) |> 
    select(-var_eta)
    # slice_max(var_eta, n=B, with_ties=FALSE) |> 
    # select(-var_eta)

d_simple_var <- tibble(date=rep(NA, B), site=rep(NA, B))
d_simple_var[1,] <- pp_grid_sort[1,]

dsize <- 1
check_idx <- 2
while (dsize <= B) {
    d_next <- pp_grid_sort[check_idx,]
    check_idx <- check_idx+1
    if (!(d_next$site %in% unique(d_simple_var$site))) {
        dsize <- dsize+1
        d_simple_var[dsize,] <- d_next
    }
}

u_simple_var <- loc_utility(
    d_simple_var, 
    pp_grid, 
    pp_parks, 
    n=20, m=1500
)

# "Retrospective" strategy:-------------------------------------------------------
# Add points with the highest variance, based on predictions from what------------
# has been seen so far. For |d|>1, use expected variance over y | d_curr----------

max_var <- function(pp_df) {
    pp_df |> 
        group_by(date, site) |> 
        summarise(var_eta=mean(var_eta), .groups="drop") |> 
        slice_max(var_eta, n=1, with_ties=FALSE) |> 
        select(-var_eta)
}

find_next_max_var <- function(d_curr, grid, parks, n) {
    nl_df <- semi_join(grid, d_curr)
    pred_loc <- anti_join(grid, d_curr)
    rep_var_grid <- map_dfr(1:n, ~{
        nl_df$pres <- rbinom(nrow(nl_df), rep(1, nrow(nl_df)), nl_df$mean_pres) # sample new y | new_loc
        mod_df <- pred_loc |>
            bind_rows(nl_df, parks)
        ft <- fit_model(f_no_corr, mod_df, response="pres", fx_prec=0.05)
        mod_df |> 
            mutate(var_eta=ft$summary.linear.predictor$sd^2) |> # get the new posterior variance
            filter(is.na(pres))
    })
    return(max_var(rep_var_grid))
}

B <- 10
d_retro0 <- max_var(pp_grid)
d_retro <- d_retro0

for (i in seq_len(B-1)) {
    d_next <- find_next_max_var(d_retro, pp_grid, pp_parks, 10)
    d_retro <- bind_rows(d_retro, d_next)
}

u_retro <- loc_utility(
    d_retro, 
    pp_grid, 
    pp_parks, 
    n=20, m=1500
)


# Prospective strategy:-----------------------------------------------------------

init_utils <- function(pred_loc, known_loc, new_loc=tibble(), n_init=10, n=10, m=5) {
    eval_pts <- pred_loc |> 
        slice_sample(n=n_init) |>
        select(X, Y, t)
        
    map_dfr(
        1:n_init,
        ~loc_utility(eval_pts[.x,], pred_loc, known_loc, new_loc, n=n, m=m, by=c("X", "Y", "t"))
        # .options=furrr_options(seed=TRUE)
    )
}

iter_bayes_opt <- function(
        f, init_utils, pred_loc, 
        max_iter=10, max_same_iter=max_iter, noise_util=0.1,
        control_km=list(trace=FALSE, pop.size=100, max.generations=50, wait.generations=10)
) {
    best_f <- min(init_utils$utility)
    nsame <- 0
    ures <- bind_rows(init_utils, tibble(utility=rep(NA, max_iter)))
    
    eval_pts <- pred_loc |>
        distinct(X, Y, t) |> 
        relocate(X, Y, t)
    
    for (i in seq_len(max_iter)) { # choose a new point to evaluate U
        ucur <- ures |> filter(!is.na(utility))
        krig <- km(
            ~1, select(ucur, X, Y, t), ucur$utility, scaling=TRUE, covtype="exp",
            noise.var=rep(noise_util, nrow(ucur)), lower=rep(0.1, 3), upper=rep(1, 3),
            control=control_km
        )
        # res <- max_AKG(
        #     krig, new.noise.var=noise_util, type="UK", lower=rep(0, 3), upper=rep(1, 3), 
        #     control=control_akg
        # )
        sub_eval_pts <- slice_sample(eval_pts, n=2500)
        eqi_pts <- map_dbl(
            1:nrow(sub_eval_pts), 
            ~EQI(c(sub_eval_pts[.x,]), krig, new.noise.var=noise_util, type="UK")
        )
        print(summary(eqi_pts))
        next_pt <- sub_eval_pts[which.max(eqi_pts),]
        next_u_eval <- f(next_pt)
        print(next_u_eval)
        if (next_u_eval$utility <= best_f) {
            best_f <- next_u_eval$utility
            nsame <- 0
        } else {
            nsame <- nsame + 1
        }
        ures[nrow(init_utils)+i,] <- next_u_eval
        if (nsame > max_same_iter)
            break
    }
    ures <- filter(ures, !is.na(utility))
    
    return(list(
        best=slice_min(ures, utility, n=1, with_ties=TRUE),
        points=ures
    ))
}

bayes_opt_design <- function() {
    
}

# first design point. Allow high levels of noise since surface is probably almost constant
utils0 <- init_utils(pp_grid, pp_parks, n_init=40, n=10, m=NULL)
f <- function(x) loc_utility(x, pp_grid, pp_parks, n=10, m=200, by=c("X", "Y", "t"))
bo <- bayes_opt(f, utils0, pp_grid, max_iter=20, max_same_iter=6)

ggplot(bo$points, aes(X, Y, size=utility)) +
    geom_point() +
    geom_point(data=u_init, col="green", shape=1) +
    facet_wrap(~t)

d_bayes_opt <- select(bo$best, -utility)
pred_loc <- anti_join(pp_grid, bo$best)
for (i in 1:1) {
    u_init <- init_utils(pred_loc, pp_parks, new_loc=d_bayes_opt, n_init=40, n=10, m=100)
    f <- function(x) loc_utility(x, pred_loc, pp_parks, d_bayes_opt, n=10, m=200, by=c("X", "Y", "t"))
    bo <- bayes_opt(f, u_init, pred_loc, max_iter=20, max_same_iter=6)
    d_bayes_opt <- bind_rows(d_bayes_opt, select(bo$best, -utility))
    pred_loc <- anti_join(pred_loc, bo$best)
}



f <- function(x) loc_utility(x, pred_loc, pp_parks, new_loc=new_loc, n=10, m=100)
bo1 <- bayes_opt(f, utils1, max_iter=15)


krig <- km(
    ~1, select(utils0, X, Y, t), utils0$utility, scaling=TRUE,
    noise.var=rep(0.1, nrow(utils0)), lower=rep(0.1, 3), upper=rep(1, 3)
)

fun <- function(x) loc_utility(pp_grid, pp_parks, x, n=10, m=100)$utility

fun(c(0.1, 0.1, 0.2))

res <- noisy.optimizer(
    optim.crit="AKG", model=krig, funnoise=fun, 
    n.ite=10, noise.var=0.1, 
    lower=rep(0, 3), upper=rep(1, 3), 
    NoiseReEstimate=FALSE, CovReEstimate=TRUE
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
krig_pred <- predict(bo1$krig, dom, type="UK", bias.correct=TRUE)

mutate(dom, kmean=krig_pred$mean) |> 
    complete(X, Y, t)


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
