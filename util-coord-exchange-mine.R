library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)
library(nngeo)

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

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

all_loc_sp <- read_sf("data-proc/parks-design-space.shp") |> 
    distinct(site, geometry)

nn_list <- nearest_neighbors(all_loc_sp, 4)

n_est <- 45
n <- 45

util_fun <- function(d) {
    utility(d, obs_mod, n=n_est, by=c("month", "site"), pred_df=pred_mod, u_only=TRUE, util_fun=util_risk_sd, risk_df=risk_mod)
}

# d_so_far <- readRDS("util-results/coord-ex-2.rds")$design[[1]]
d_so_far <- ce_res2$d
ce_res3 <- coordinate_exchange(util_fun, pred_mod, 5, nn_list, max_iter=100, d_so_far=d_so_far, verbose=TRUE)

res <- tibble_row(risk_sd=ce_res$u_approx, design=list(ce_res$d), num_loc=nrow(ce_res$d))
saveRDS(res, "util-results/risk-sd/coord-ex-4.rds")

res <- map_dfr(1:4, ~{
    res_d <- readRDS(paste0("util-results/risk-sd/coord-ex-", .x, ".rds"))
    res_d
})

saveRDS(res, "util-results/risk-sd/util-coord-ex.rds")
