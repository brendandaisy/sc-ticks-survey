library(INLA)
library(purrr)
library(tidyverse)
library(sf)

library(tidyverse)
library(tidygraph)
library(ggraph)
library(furrr)

inla.setOption(inla.mode="classic", num.threads="16:1")

source("other-helpers.R")
source("utility-helpers.R")

parks_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

all_data <- append_pred_grid(parks_data) |> st_drop_geometry()

grid_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
parks_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

num_d_rep <- 3 # number of designs, for each design size
num_loc <- c(1, 10, 20)
ns <- c(10, 20, 40)
set.seed(155)
plan(future::multicore(workers=8))

# a random, first design of 3 points to test with

design_err <- function(num_loc, ns, n_reps=1, sel=sel_list_inla(parks_mod)) {
    d <- grid_mod |> 
        select(date, site) |> 
        slice_sample(n=num_loc)
    
    d_df <- semi_join(grid_mod, d, by=c("date", "site"))
    new_df <- prep_new_data(parks_mod, d_df, scale=FALSE)
    
    # get one sample of [ynew | y, d]
    pred_idxs <- which(is.na(new_df$pres))
    pred_mat <- rpost_predict(new_df, pred_idxs, n=1)
    new_df$pres[pred_idxs] <- rbinom(nrow(d_df), rep(1, nrow(d_df)), pred_mat[,1])
    # error for a particular ynew
    urep_yfix <- map_dfr(1:n_reps, ~tibble_row(!!!util_rep(new_df, sel)))
    
    # error for proper utility (many samples of ynew)
    ns <- rep(ns, each=n_reps) # for each n, do n_reps
    urep_yrep <- map_dfr(ns, ~mutate(tibble_row(!!!utility(d, parks_mod, .x, grid_mod, u_only=TRUE)), n=.x))
    
    tibble_row(
        design=list(d),
        urep_yfix=list(urep_yfix), 
        urep_yrep=list(urep_yrep)
    )
}

d_sizes <- rep(num_loc, each=num_d_rep)
future_iwalk(
    d_sizes, 
    ~saveRDS(design_err(.x, ns=ns, n_reps=10), paste0("mcmc-err", .y, ".rds")),
    .options=furrr_options(seed=TRUE),
    .progress=TRUE
)

res <- map_dfr(seq_along(d_sizes), ~{
    res_d <- readRDS(paste0("mcmc-err", .x, ".rds"))
    # file.remove(paste0("mcmc-err", .x, ".rds"))
    res_d
})

res$urep_yfix

saveRDS(res, "mcmc-err-results.rds")

res$urep_yrep

err_results <- res |> 
    mutate(d_id=1:n()) |> 
    unnest(urep_yrep) |> 
    pivot_longer(Dfixed:Eres, "criteria", values_to="utility") |> 
    group_by(n, d_id, criteria) |> 
    summarise(sd=sd(utility), snr=abs(mean(utility) / sd))

ggplot(err_results, aes(as.character(n), sd)) +
    geom_point() +
    geom_boxplot() +
    facet_wrap(~criteria, scales="free") +
    labs(x="N", y="Signal to Noise")
