library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

nearest_neighbors <- function(loc_sp, k) {
    nn_list <- loc_sp |> 
        st_nn(loc_sp, k=k+1) |> 
        map(\(nn) loc_sp$site[nn[-1]])
    
    names(nn_list) <- loc_sp$site
    return(nn_list)
}

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
            for (m in (m_curr+1):12) { # cycle through dates until utility is not improved
                d_prop$month[i] <- if ((m_curr + m) == 12) 12 else (m_curr + m) %% 12
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
                print(d_prop)
                print(u_prop)
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

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-dopt-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
# risk_mod <- df_list$risk_grid

n_est <- 10
n <- 50

util_fun <- function(d) {
    utility(d, obs_mod, n=30, by=c("month", "site"), pred_df=pred_mod, u_only=TRUE, sel=sel_list_inla(obs_mod))
}

d_so_far <- tibble()
nn_list <- nearest_neighbors(all_loc_sp, 2)
coordinate_exchange(util_fun, pred_mod, 5, nn_list, max_iter=25, verbose=FALSE)
for (i in 1:2) {
    sa_res <- simulated_annealing(pred_mod, obs_mod, 5, iter=100, alpha=1.3, T0=0.02, n_est=n_est, d_so_far=d_so_far, info_iter=5, restart=50)
    d_so_far <- sa_res$d
    u_res <- utility(sa_res$d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)
    res <- tibble_row(!!!u_res, num_loc=nrow(d_so_far))
    saveRDS(res, paste0("util-results/risk-sd/sim-ann3-", i, ".rds"))
}
# d_init <- res$design
# for (i in 3) {
#     sa_res <- simulated_annealing(pred_mod, obs_mod, 15, iter=20, alpha=1.7, T0=0.01, n_est=n_est, d_init=d_init[[i]], info_iter=5)
#     d_so_far <- sa_res$d
#     u_res <- utility(sa_res$d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)
#     res <- tibble_row(!!!u_res, num_loc=nrow(d_so_far))
#     saveRDS(res, paste0("util-results/util-sim-ann-alpha1pt5-T0pt2-mod-", i, ".rds"))
# }

res <- map_dfr(1:2, ~{
    res_d <- readRDS(paste0("util-results/risk-sd/sim-ann2-", .x, ".rds"))
    res_d
})

saveRDS(res, "util-results/util-sim-ann.rds")