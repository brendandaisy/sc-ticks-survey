# --------------------------------------------------------------------------------
# util-random-sampling.R----------------------------------------------------------
# script to sample designs uniformly at random, eval their utility, and save------
# the results---------------------------------------------------------------------
# --------------------------------------------------------------------------------
library(INLA)
library(tidyverse)
library(sf)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

# df_list <- readRDS("data-proc/bo-dopt-dfs.rds")
# obs_mod <- df_list$obs
# pred_mod <- df_list$pred

num_loc <- seq(5, 20, 5)
reps <- 3
expg <- expand_grid(num_loc=num_loc, rep=1:reps)
n <- 50

# Random selection, baseline design strategy--------------------------------------
u_random <- function(num_loc, i, app_ind=0) {
    d <- pred_mod |> 
        select(date, site) |> 
        slice_sample(n=num_loc)
    
    res <- tibble_row(
        !!!utility(d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod),
        num_loc=num_loc
    )
    print(res)
    saveRDS(res, paste0("util-results/risk-sd/util-random/util-random", i+app_ind, ".rds"))
}

iwalk(expg$num_loc, u_random, app_ind=81)

res <- map_dfr(1:93, ~readRDS(paste0("util-results/risk-sd/util-random/util-random", .x, ".rds")))

saveRDS(res, "util-results/risk-sd/util-random.rds")
