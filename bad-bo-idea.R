library(tidyverse)
library(INLA)
library(sf)

source("other-helpers.R")

poly_grid <- read_sf("geo-files/canopy-poly-4km.shp") |> 
    rename(tree_canopy=tr_cnpy) |> 
    mutate(site=1:n())

pred_model_data <- map_dfr(1:10, ~mutate(poly_grid, month=.x)) |> 
    mutate(y=1.3+tree_canopy*1.2+-.5*month+rnorm(n(), sd=2))

nobs <- nrow(pred_model_data)

pred_model_data$y[sample(nobs, nobs-10)] <- NA

grid_nb <- poly2nb(poly_grid, row.names=poly_grid$site)
grid_mat <- as(nb2mat(grid_nb, style="B"), "Matrix")

plot(pd_nb, st_coordinates(poly_grid))

loc <- st_coordinates(pred_model_data)

fit <- inla(
    y ~ 1 + f(site, model="besag", graph=grid_mat, group=month, control.group=list(model="rw1")),
    data=pred_model_data,
    control.predictor=list(compute=TRUE)
)
