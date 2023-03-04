library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)

source("other-helpers.R")
source("utility-helpers.R")
dir.create("util-results", showWarnings=FALSE)

# set.seed(203)
# inla.setOption(inla.mode="classic", num.threads="16:1")
# plan(future::multicore(workers=8))

init_d_cand <- function(pred_df) {
    pred_df |> 
        group_by(date, site) |> 
        summarise(var_eta=mean(var_eta), .groups="drop")
}

local_utility <- function(d, known_df, pred_df=NULL, copy=FALSE) {
    # TODO: revisit at some point whether known_df can be removed, or can be used for only one tick species
    # probably can, but confusing especially with the dummy matrix
    d_df <- if (!is.null(pred_df)) inner_join(pred_df, d, by=c("date", "site"), copy=copy) else d
    d_df <- bind_rows(d_df, known_df)
    
    X <- dummyVars(
        pres ~ -1 + tick_class + tick_class:tree_canopy + tick_class:elevation +
            #^ can't include land cover because different designs have diff levels,
            # and penalizes designs with new levels incorrectly, but also can't just keep all levels since singular
            tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh,
        data=d_df,
        fullRank=TRUE
    )
    X <- predict(X, d_df)
    
    I <- t(X) %*% diag(d_df$mean_pres * (1 - d_df$mean_pres)) %*% X
    return(log(det(I))) # want to MAXIMIZE information to minimize uncertainty
}

sim_ann_local <- function(
        pred_df, known_df, num_loc=1, iter=1, alpha=1, T0=1, weight_by=NULL
) {
    # assign each date, site a weight based on avg variance of each species
    # this is the set new designs points will be chosen from:
    d_cand <- init_d_cand(pred_df) 
    d_curr <- d_best <- slice_sample(d_cand, n=num_loc, weight_by={{weight_by}})

    u_curr <- u_best <- local_utility(d_curr, known_df, pred_df)
    i_best <- 0
    
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    for (i in 1:iter) {
        # choose a new seach point. Reusing a point in d_curr IS allowed
        next_pt <- slice_sample(d_cand, n=1, weight_by={{weight_by}})
        
        # delete a random point and add the new one
        if (num_loc > 1) {
            d_prop <- slice_sample(d_curr, n=nrow(d_curr)-1) |>
                bind_rows(next_pt)
        }
        else
            d_prop <- next_pt
        
        u_prop <- local_utility(d_prop, known_df, pred_df)
        if (u_prop > u_best) { # if better than global best, update that info
            i_best <- i
            d_curr <- d_best <- d_prop
            u_curr <- u_best <- u_prop
        } else if (u_prop > u_curr) { # if better than current design, accept but don't update global best
            d_curr <- d_prop
            u_curr <- u_prop
        } else if (runif(1) < exp((log10(u_prop) - log10(u_curr))/T_sched[i])) { # or accept anyway sometimes
            d_curr <- d_prop
            u_curr <- u_prop
        }
    }
    return(list(d=d_best, u_local=u_best, at=i_best))
}

# Data preparation----------------------------------------------------------------
parks_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |>
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

fit_all <- readRDS("inla-fit-all-data.rds")

all_data <- append_pred_grid(parks_data) |>
    st_drop_geometry() |> 
    mutate(
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd,
        var_eta=fit_all$summary.linear.predictor$sd^2
    )

grid_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
parks_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

n <- 50
num_loc <- c(1, seq(5, 20, 5))

res <- map_dfr(num_loc, ~{
    sa_res <- sim_ann_local(grid_mod, parks_mod, .x, iter=1000, alpha=1.5, T0=0.75)
    u_res <- utility(sa_res$d, parks_mod, n=n, full_df=grid_mod)
    print(u_res$dopt)
    tibble_row(!!!u_res, dopt_local=sa_res$u_local, num_loc=.x)
})
