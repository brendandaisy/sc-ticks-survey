library(INLA)
library(purrr)
library(tidyverse)
library(sf)
library(furrr)

inla.setOption(inla.mode="classic", num.threads="16:1")

source("other-helpers.R")
source("utility-helpers.R")

dir.create("util-results", showWarnings=FALSE)

set.seed(202)
plan(future::multicore(workers=8))

# Data preparation----------------------------------------------------------------
parks_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

all_data <- append_pred_grid(parks_data) |> 
    st_drop_geometry()

grid_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
parks_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

num_loc <- c(1, seq(5, 20, 5))
reps <- 15
expg <- expand_grid(num_loc=num_loc, rep=1:reps)
n <- 50

# Random selection, baseline design strategy--------------------------------------
u_random <- function(num_loc, i) {
    d <- grid_mod |> 
        select(date, site) |> 
        slice_sample(n=num_loc)
    
    res <- tibble_row(!!!utility(d, parks_mod, n=n, full_df=grid_mod), num_loc=num_loc)
    saveRDS(res, paste0("util-results/util-random", i, ".rds"))
}

future_iwalk(expg$num_loc, u_random, .options=furrr_options(seed=TRUE), .progress=TRUE)

res <- map_dfr(1:nrow(expg), ~{
    res_d <- readRDS(paste0("util-results/util-random", .x, ".rds"))
    # file.remove(paste0("mcmc-err", .x, ".rds"))
    res_d
})

saveRDS(res, "util-results/util-random.rds")
