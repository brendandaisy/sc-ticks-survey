# --------------------------------------------------------------------------------
# util-one-offs.R-----------------------------------------------------------------
# evaluate utility of the one off/repeated sampling designs-----------------------
# --------------------------------------------------------------------------------
library(INLA)
library(purrr)
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

n <- 50

# Visit all 30 parks again in December/June
d_parks_dec <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-12-01"))

u_parks_dec <- utility(d_parks_dec, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)

d_parks_jun <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-06-01"))

u_parks_jun <- utility(d_parks_jun, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)

# Reuse the same schedule from 2021
d_resamp <- filter(obs_mod, lubridate::year(date) == 2021)
d_resamp$pres <- NA
u_resamp <- utility(d_resamp, obs_mod, n=n, util_fun=util_risk_sd, risk_df=risk_mod)

u_one_offs <- tibble(
    strat=c("    Revisit in\n    June (30)", "    Revisit in\n    December (30)", "    Repeat 2021\n    schedule (111)"), 
    design=list(distinct(d_parks_jun, date, site), distinct(d_parks_dec, date, site), distinct(d_resamp, date, site)),
    utility=c(u_parks_jun$risk_sd, u_parks_dec$risk_sd, u_resamp$risk_sd)
)

saveRDS(u_one_offs, "util-results/risk-sd/util-one-offs.rds")