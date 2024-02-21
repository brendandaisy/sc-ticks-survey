library(tidyverse)
library(sf)

parks_obs <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE)

## First save the 2 park dfs for D-optimality criterion separately
## This is mainly to not penalize for land_cover classes that aren't in the parks data, but also 
## is rescaled with different constants (since no grid data)
all_data <- append_pred_data(parks_obs) |> 
    st_drop_geometry()

pred_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
obs_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

saveRDS(list(obs=obs_mod, pred=pred_mod), "data-proc/bo-dopt-dfs.rds")

##
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

# only include visits to sites with risk > 0.7 for 3+ months and/or species per year
risk_grid_mod <- all_data |> 
    filter(is.na(pres) & str_detect(site, "g\\d+")) |> 
    group_by(date, site) |> 
    mutate(at_risk=sum(mean_pres > 0.75)) |> 
    ungroup() |> 
    filter(at_risk >= 1)

distinct(risk_grid_mod, site)

saveRDS(list(obs=obs_mod, pred=pred_mod, risk_grid=risk_grid_mod), "data-proc/bo-risk-sd-dfs.rds")