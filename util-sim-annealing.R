# --------------------------------------------------------------------------------
# util-sim-annealing.R------------------------------------------------------------
# search for an optimal design using a Simulated Annealing algorithm, and---------
# save the results----------------------------------------------------------------
# --------------------------------------------------------------------------------
library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

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
#' @param info_iter print progress after this many iterations has passed
#' @param restart optionally, an iteration to reset back to the best utility found so far
#'
#' @return a `list` contain the best design, its approx utility, and the iteration it was found in
simulated_annealing <- function(
        pred_df, known_df, num_loc=1, iter=1, alpha=1, T0=1, n_est=1, 
        d_so_far=tibble(), d_init=NULL, weight_by=var_eta, info_iter=10, restart=NULL
) {
    # assign each date, site a weight based on avg variance of each species
    # this is the set new designs points will be chosen from:
    d_cand <- init_d_cand(pred_df)
    
    if (is.null(d_init)) # then chose initial design randomly
        d_curr <- d_best <- slice_sample(d_cand, n=num_loc, weight_by={{weight_by}})
    else
        d_curr <- d_best <- d_init
    
    u_curr <- u_best <- utility(
        bind_rows(d_curr, d_so_far), known_df, 
        n=n_est, pred_df=pred_df, util_fun=util_risk_sd, risk_df=risk_mod, u_only=TRUE
    )
    i_best <- 0
    
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    for (i in 1:iter) {
        # choose a new seach point. Reusing a point in d_curr IS allowed
        next_pt <- slice_sample(d_cand, n=1, weight_by={{weight_by}})
        
        # delete a random point and add the new one
        d_prop <- slice_sample(d_curr, n=max(nrow(d_curr)-1, 1)) |>
            bind_rows(next_pt)
        
        u_prop <- utility(
            bind_rows(d_prop, d_so_far), known_df, 
            n=n_est, pred_df=pred_df, util_fun=util_risk_sd, risk_df=risk_mod, u_only=TRUE
        )
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
        if (i %% info_iter == 0) # print update after iter
            print(paste0("Current U=", round(u_curr, 3), " at i=", i, "/", iter))
        if (!is.null(restart) && i == restart) {
            d_curr <- d_best
            u_curr <- u_best
        }
    }
    return(list(d=bind_rows(d_best, d_so_far), u_approx=u_best, at=i_best))
}

init_d_cand <- function(pred_df) {
    pred_df |> 
        group_by(date, site) |> 
        summarise(var_eta=mean(var_eta), .groups="drop")
}

# helper to plot tuning of alpha and T0 for SA
sa_acceptance_prob <- function(u_prop, u_curr, alpha, T0, iter) {
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    probs <- exp((log10(u_prop) - log10(u_curr)) / T_sched)
    
    ggplot(enframe(probs), aes(name, value)) +
        geom_point() +
        labs(x="iteration", y="prob.")
}

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

n_est <- 50 # can make this lower for faster search
n <- 50

d_so_far <- tibble()
for (i in 1:4) {
    sa_res <- simulated_annealing(pred_mod, obs_mod, 5, iter=150, alpha=1.3, T0=0.02, n_est=n_est, d_so_far=d_so_far, info_iter=5, restart=50)
    d_so_far <- sa_res$d
    u_res <- utility(sa_res$d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)
    res <- tibble_row(!!!u_res, num_loc=nrow(d_so_far))
    saveRDS(res, paste0("util-results/risk-sd/sim-ann3-", i, ".rds"))
}

# combine the results for each number of visits to a single file
res <- map_dfr(1:4, ~{
    res_d <- readRDS(paste0("util-results/risk-sd/sim-ann-", .x, ".rds"))
    res_d
})
 
saveRDS(res, "util-results/risk-sd/util-sim-ann.rds")
