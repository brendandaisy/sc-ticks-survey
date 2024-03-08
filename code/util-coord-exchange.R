# --------------------------------------------------------------------------------
# util-coord-exchange.R-----------------------------------------------------------
# search for an optimal design using an exchange algorithm and save the results---
# --------------------------------------------------------------------------------
library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)
library(nngeo)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

#' coordinate_exchange
#'
#' @param util_fun a wrapper around `utility` accepting a single argument for the design
#' Allows fixing utility settings such as `n` beforehand 
#' @param pred_df the prediction data, or where designs can choose to visit. Requires all columns present
#' @param num_loc the number of locations to add in this round of the design
#' @param nn_list an adjacency list of names as candidate locations and values as locations considered neighbors
#' @param max_iter maximum number of allowed utility evaluations. Note convergence of the algorithm can lead
#' to fewer iterations in total
#' @param d_so_far optionally, a previously found design, to add `num_loc` more on to. These points will NOT be changed
#' @param d_init optionally, a pre specified design to start searching from. These points WILL be changed
#' @param verbose prints utility so far after each "sweep"
#'
#' @return a `list` contain the best design and its approx utility 
coordinate_exchange <- function(
        util_fun, pred_df, num_loc, nn_list, max_iter=10,
        d_so_far=tibble(), d_init=NULL, verbose=FALSE
) {
    if (is.null(d_init)) # then chose initial design randomly
        d_curr <- pred_df |> distinct(month, site) |> slice_sample(n=num_loc)
    else
        d_curr <- d_init
    
    u_curr <- util_fun(bind_rows(d_curr, d_so_far))
    iter <- 1
    while (iter <= max_iter) {
        changed <- FALSE
        for (i in 1:num_loc) { # for each visit in current design
            d_prop <- d_curr # reset d_prop if last site from visit i-1 rejected
            m_curr <- d_curr$month[i]  
            for (m in 1:11) { # cycle through dates until utility is not improved
                d_prop$month[i] <- if ((m_curr + m) == 12) 12 else (m_curr + m) %% 12
                u_prop <- util_fun(bind_rows(d_prop, d_so_far))
                iter <- iter + 1
                if (u_prop < u_curr)
                    break
                d_curr <- d_prop
                u_curr <- u_prop
                changed <- TRUE
            }
            for (m in 1:11) { # cycle backwards through dates until utility is not improved
                d_prop$month[i] <- if ((m_curr - m) == 0) 12 else (m_curr - m) %% 12
                u_prop <- util_fun(bind_rows(d_prop, d_so_far))
                iter <- iter + 1
                if (u_prop < u_curr)
                    break
                d_curr <- d_prop
                u_curr <- u_prop
                changed <- TRUE
            }
            s_curr <- d_curr$site[i]
            d_prop <- d_curr # reset d_prop if the last month was rejected
            for (nn in nn_list[[s_curr]]) {
                d_prop$site[i] <- nn
                u_prop <- util_fun(bind_rows(d_prop, d_so_far))
                iter <- iter + 1
                if (u_prop > u_curr) {
                    d_curr <- d_prop
                    u_curr <- u_prop
                    changed <- TRUE
                }
            }
        }
        if (verbose) {
            print("Current d:")
            print(d_curr)
            print(paste0("Current U: ", round(u_curr, 3), " at i=", iter, "/", max_iter))
        }
        if (!changed)
            break
    }
    return(list(d=bind_rows(d_curr, d_so_far), u_approx=u_curr))
}

# find a nearest neighbors adjacency list using `st_nn`
# k is the number of neighbors (self-loops are removed)
nearest_neighbors <- function(loc_sp, k) {
    nn_list <- loc_sp |> 
        st_nn(loc_sp, k=k+1) |> 
        map(\(nn) loc_sp$site[nn[-1]])
    
    names(nn_list) <- loc_sp$site
    return(nn_list)
}

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds") # change to bo-dopt-dfs for the first criterion
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

all_loc_sp <- read_sf("data-proc/parks-design-space.shp") |> 
    distinct(site, geometry)

nn_list <- nearest_neighbors(all_loc_sp, 4)

n_est <- 50
n <- 50

util_fun <- function(d) { # would also change for the first criterion
    utility(d, obs_mod, n=n_est, by=c("month", "site"), pred_df=pred_mod, u_only=TRUE, util_fun=util_risk_sd, risk_df=risk_mod)
}

res <- vector("list", 4)
d_so_far <- tibble()
for (i in 1:4) {
    ex_res <- coordinate_exchange(util_fun, pred_mod, 5, nn_list, max_iter=150, d_so_far=d_so_far, verbose=TRUE)
    d_so_far <- ex_res$d
    u_res <- utility(ex_res$d, obs_mod, n=n, by=c("month", "site"), pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)
    res[[i]] <- list(u_res=u_res, d=d_so_far)
}
