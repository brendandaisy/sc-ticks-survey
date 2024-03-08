# --------------------------------------------------------------------------------
# prep-bedopt-dfs.R---------------------------------------------------------------
# A convenience script to produce and save data that will be used during the BED--
# analysis. To get utlity of a proposed design, load the files matching the-------
# utility criterion. Files include initial collection data with appropriate-------
# scaling of environmental variables, future possible visits, and the high risk---
# grid locations (for the second criterion)---------------------------------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(sf)

parks_obs <- read_parks_sf()

# For the D-optimality criterion, save the 2 park dfs of initial visits and future visits
# This is mainly to not penalize for land_cover classes that aren't in the parks data, but also 
# covariates are rescaled with different values (since no grid data)
all_data <- append_pred_data(parks_obs, drop_new_lcc=TRUE) |> 
    st_drop_geometry()

pred_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
obs_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

saveRDS(list(obs=obs_mod, pred=pred_mod), "data-proc/bo-dopt-dfs.rds")

# For the high risk location criterion, save 3 dfs. The initial and future visits with different scaling, and 
# the high risk grid locations used to calc utility
all_data <- append_pred_data(parks_obs, pred_grid=prep_pred_grid(parks_obs)) |> 
    st_drop_geometry()

fit_all <- fit_model(formula_bed(), all_data, fx_prec=0.2)

all_data <- mutate(
    all_data,
    mean_pres=fit_all$summary.fitted.values$mean,
    sd_pres=fit_all$summary.fitted.values$sd,
    var_eta=fit_all$summary.linear.predictor$sd^2
)

pred_mod <- filter(all_data, is.na(pres) & !str_detect(site, "g\\d+")) # final not yet observed data for fitting models
obs_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

# only include visits to sites with risk > 0.75 for at least one species
risk_grid_mod <- all_data |> 
    filter(is.na(pres) & str_detect(site, "g\\d+")) |> 
    group_by(date, site) |> 
    mutate(at_risk=sum(mean_pres > 0.75)) |> 
    ungroup() |> 
    filter(at_risk >= 1)

saveRDS(list(obs=obs_mod, pred=pred_mod, risk_grid=risk_grid_mod), "data-proc/bo-risk-sd-dfs.rds")