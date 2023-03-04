library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)

source("other-helpers.R")
source("utility-helpers.R")
dir.create("util-results", showWarnings=FALSE)

set.seed(203)
inla.setOption(inla.mode="classic", num.threads="16:1")
plan(future::multicore(workers=8))

init_d_cand <- function(pred_df) {
    pred_df |> 
        group_by(date, site) |> 
        summarise(var_eta=mean(var_eta), .groups="drop")
}

#' simulated_annealing
#'
#' @param pred_df the prediction data, or where designs can choose to observe. Requires all columns present
#' @param known_df data of existing observations, that cannot be part of new design but are added to each design
#' @param num_loc the number of locations to add in this round of the design
#' @param iter the number of steps to run the SA algorithm
#' @param alpha curvature of cooling schedule. `alpha` > 1 means faster cooling (less acceptance of worse moves)
#' @param T0 cooling magnitude. The initial order of mag. worse a prop can be, for an acceptance prob of exp(-1). `T0`>1 means more exploring
#' @param n_est number of MC iterations to estimate utility
#' @param d_so_far optionally, a previously found design, to add `num_loc` more on to. These points will NOT be changed
#' @param d_init optionally, a pre specified design to start searching from. These points WILL be changed
#' @param weight_by a variable to weight by when choosing new design points
#'
#' @return a `list` contain the best design, its approx utility, and the iteration it was found in
simulated_annealing <- function(
        pred_df, known_df, num_loc=1, iter=1, alpha=1, T0=1, n_est=1, 
        d_so_far=tibble(), d_init=NULL, weight_by=var_eta
) {
    # assign each date, site a weight based on avg variance of each species
    # this is the set new designs points will be chosen from:
    d_cand <- init_d_cand(pred_df) 
    
    if (is.null(d_init)) # then chose initial design randomly
        d_curr <- d_best <- slice_sample(d_cand, n=num_loc, weight_by={{weight_by}})
    else
        d_curr <- d_best <- d_init
    u_curr <- u_best <- utility(bind_rows(d_curr, d_so_far), known_df, n=n_est, full_df=pred_df, u_only=TRUE)
    i_best <- 0
    
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    for (i in 1:iter) {
        # choose a new seach point. Reusing a point in d_curr IS allowed
        next_pt <- slice_sample(d_cand, n=1, weight_by={{weight_by}})
        
        # delete a random point and add the new one
        d_prop <- slice_sample(d_curr, n=max(nrow(d_curr)-1, 1)) |>
            bind_rows(next_pt)
        
        u_prop <- utility(bind_rows(d_prop, d_so_far), known_df, n=n_est, full_df=pred_df, u_only=TRUE)
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
        if (i %% 10 == 0) # print update after 10 iter
            print(paste0("Current U=", round(u_curr, 3), " at i=", i, "/", iter))
    }
    return(list(d=bind_rows(d_best, d_so_far), u_approx=u_best, at=i_best))
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

n_est <- 30
n <- 50

d_so_far <- tibble()
for (i in 1:4) {
    sa_res <- simulated_annealing(grid_mod, parks_mod, 5, iter=100, alpha=1.5, T0=0.1, n_est=n_est, d_so_far=d_so_far)
    d_so_far <- sa_res$d
    u_res <- utility(sa_res$d, parks_mod, n=n, full_df=grid_mod)
    res <- tibble_row(!!!u_res, num_loc=nrow(d_so_far))
    saveRDS(res, paste0("util-results/util-sim-ann-alpha1pt5-T0pt1", i, ".rds"))
}

# res <- map_dfr(1:4, ~{
#     res_d <- readRDS(paste0("util-results/util-sim-ann", .x, ".rds"))
#     # file.remove(paste0("mcmc-err", .x, ".rds"))
#     res_d
# })
# 
# saveRDS(res, "util-results/util-sim-ann.rds")
