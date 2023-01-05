library(INLA)
library(tidyverse)
library(sf)

source("other-helpers.R")
source("utility-helpers.R")

parks_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

all_data <- append_pred_grid(parks_data) |> st_drop_geometry()
fit_all <- readRDS("inla-fit-all-data.rds")

post_pred_all <- all_data |>
    mutate(
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd,
        var_eta=fit_all$summary.linear.predictor$sd^2
    ) |> 
    st_drop_geometry() # drop here, since need to at some point and don't need to project to R^3

# pp_grid <- filter(post_pred_all, is.na(pres))
# pp_parks <- filter(post_pred_all, !is.na(pres))
grid_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
parks_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

# a random, first design of 3 points to test with
set.seed(125)

d3 <- grid_mod |> 
    select(date, site) |> 
    slice_sample(n=3)

d3_df <- semi_join(grid_mod, d3)
new_df <- prep_new_data(parks_mod, d3_df, scale=FALSE)
pred_idxs <- which(is.na(new_df$pres))

pred_mat <- rpost_predict(new_df, pred_idxs, n=1)
new_df$pres[pred_idxs] <- rbinom(nrow(d3_df), rep(1, nrow(d3_df)), pred_mat[,1])

urep_yfix <- map_dbl(1:100, ~util_bd_rep(new_df, sel_list_inla(parks_mod)))

mean(urep_yfix)

new_df <- prep_new_data(parks_mod, d3_df, scale=FALSE)

urep_n10 <- map_dbl(1:20, ~utility(d3, parks_mod, 10, grid_mod, u_only=TRUE))
