# --------------------------------------------------------------------------------
# other-helpers.R-----------------------------------------------------------------
# functions for reading in data, fitting INLA models, and other misc. tasks-------
# --------------------------------------------------------------------------------
library(tidyverse)
library(INLA)

# Initial data preparation--------------------------------------------------------

#' read_parks_sf
#' read in the initial collections data and remove the ESRI abbreviated names
#'
#' @param f file string
#'
#' @return an `sf` simple feature collection
read_parks_sf <- function(f="data-proc/parks-observed.shp") {
    ret <- read_sf(f) |> 
        rename(
            tick_class=tck_cls, land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, 
            max_temp=max_tmp, precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        mutate(land_cover=fct_drop(factor(land_cover, levels=land_cover_names())))
    
    mutate(ret, tcnum=case_when(
        tick_class == "Amblyomma americanum" ~ 1L, 
        tick_class == "Dermacentor variabilis" ~ 2L, 
        TRUE ~ 3L
    )) |> 
        arrange(tcnum)
}


#' rescale_covars
#'
#' @param df A tibble containing environmental covariates to rescale
rescale_covars <- function(df) {
    mutate(df, across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)))
}

#' append_pred_data
#' For the BED analysis, append the space of possible future park visits to initial
#' collections data. These visits have pres=NA so if they are chosen by a design, their
#' outcome is predicted. Optionally, `pred_grid` allows a regular grid of prediction points
#' to be added similarly (useful for the second design criteria and Figure 3)
#'
#' @param obs_df initial collections data
#' @param f file name
#' @param pred_grid output of `prep_pred_grid`
#' @param drop_new_lcc whether land cover classes that do not appear in the initial data should be dropped
#'
#' @return An sf collection with environmental variables centered and scaled
append_pred_data <- function(obs_df, f="data-proc/parks-design-space.shp", pred_grid=NULL, drop_new_lcc=TRUE) {
    pred_df <- read_sf(f) |> 
        rename(
            land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
            precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        mutate(land_cover=land_cover_labels(land_cover)) |> 
        relocate(land_cover, tree_canopy, elevation, min_temp, max_temp, precipitation, jan_min_temp, mean_rh, .after=date) |> 
        select(-min_temp)
    
    # Add mising pres for each tick species for the future visits
    pred_df <- pred_df |> 
        mutate(data=list(tibble(tick_class=unique(obs_df$tick_class), pres=NA, tcnum=1:3))) |> 
        unnest(c(data))
    
    # Optionally add the pred grid. WARNING: this will result in different scaling than if pred grid were not included
    if (!is.null(pred_grid))
        pred_df <- bind_rows(pred_df, pred_grid)
    
    if (drop_new_lcc) {
        pred_df <- filter(pred_df, land_cover %in% levels(obs_df$land_cover)) |> 
            mutate(land_cover=fct_drop(land_cover))
    }
    # all data will now be scaled, and is the reference covar values for all subdesigns
    return(rescale_covars(bind_rows(obs_df, pred_df)))
}

#' prep_pred_grid
#'
#' @param parks_data initial collection data. Really only used to get tick_class labels
#' @param f file name
#'
#' @return An sf collection matching the same format as `read_parks_sf`
prep_pred_grid <- function(parks_data, f="geo-files/covar-grid-16km.shp") {
    covar_grid <- read_sf(f) |> 
        rename(
            land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
            precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        mutate(land_cover=land_cover_labels(land_cover)) |> 
        select(-min_temp)
    
    # Add a unique site label to each grid location (i.e. const. over time)
    covar_grid |> 
        group_by(geo=as.character(geometry)) |> 
        mutate(site=str_c("g", cur_group_id())) |> 
        ungroup() |> 
        select(-geo) |> 
        mutate(data=list(tibble(tick_class=unique(parks_data$tick_class), pres=NA, tcnum=1:3))) |> 
        unnest(c(data))
}

# Helpers for working with INLA objects-------------------------------------------

#' fit_model
#'
#' @param formula The model structure specified with an INLA formula
#' @param data A tibble with data to fit to
#' @param fx_prec prior precision on the environmental effects
#' @param control_compute control options such as whether to compute DIC
#' @param ... other arguments to pass to INLA call
#'
#' @return An `inla` fit object
fit_model <- function(formula, data, fx_prec, control_compute=list(mlik=FALSE), ...) {
    if (is(data, "sf"))
        data <- st_drop_geometry(data) # convert it since INLA can be wierd with sf tibbles
    
    inla(
        formula, data=data, family="binomial",
        control.fixed = list(
            expand.factor.strategy="inla", # nec. since alternative removes unused lvls
            prec=fx_prec
        ),
        control.compute=control_compute,
        control.predictor=list(link=1, compute=TRUE), 
        ...
    )
}

# the preferred model for performing BED
formula_bed <- function(...) {
    prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.1)))

    pres ~ -1 + tick_class + land_cover + tree_canopy + elevation + max_temp + 
        precipitation + jan_min_temp + mean_rh + 
        f(month, model = "ar1", hyper = prec_pri, group = tcnum, control.group = list(model = "iid", hyper = prec_pri)) +
        f(site, model = "iid", hyper = prec_pri, group = tcnum, control.group = list(model = "iid", hyper = prec_pri))
}

# Helpers for relabelling and preparing fixed effects for analysis----------------
summ_fx <- function(fx_marg) {
    imap_dfr(fx_marg, ~{
        ci5 <- inla.hpdmarginal(0.5, .x)
        ci95 <- inla.hpdmarginal(0.95, .x)
        tibble(
            var = .y,
            m = inla.zmarginal(.x, TRUE)$mean,
            lo5 = ci5[,1], hi5 = ci5[,2],
            lo95 = ci95[,1], hi95 = ci95[,2]
        )
    })
}

fx_labels <- function(fx) {
    fx |> 
        str_remove("tick_class|land_cover") |> 
        fct_relevel("Dermacentor variabilis", "Ixodes spp.", "Amblyomma americanum", after=Inf) |> 
        fct_relevel("tree_canopy", "elevation", "jan_min_temp", "max_temp", "precipitation", "mean_rh", after=0) |> 
        fct_recode(
            "Intercept, D. variabilis"="Dermacentor variabilis", 
            "Intercept, Ixodes spp."="Ixodes spp.", 
            "Intercept, A. americanum"="Amblyomma americanum",
            "Tree canopy"="tree_canopy", 
            "Elevation"="elevation", 
            "Jan. min. temperature"="jan_min_temp", 
            "Max temperature"="max_temp", 
            "Precipitation"="precipitation", 
            "Relative humidity"="mean_rh"
        )
}

land_cover_labels <- function(lc) {
    lv <- c(11, 12, 21:24, 31, 41:43, 51, 52, 71:74, 81, 82, 90, 95) |> as.character()
    lab <- land_cover_names()
    ret <- factor(as.character(lc), levels=lv, labels=lab)
    return(fct_drop(ret))
}

land_cover_names <- function() {
    c(
        "Open Water",
        "Perennial Ice/Snow",
        "Developed, Open Space",
        "Developed, Low Intensity",
        "Developed, Medium Intensity",
        "Developed, High Intensity",
        "Barren Land",
        "Deciduous Forest",
        "Evergreen Forest",
        "Mixed Forest",
        "Dwarf Scrub",
        "Shrub/Scrub",
        "Grassland/Herbaceous",
        "Sedge/Herbaceous",
        "Lichens",
        "Moss",
        "Pasture/Hay",
        "Cultivated Crops",
        "Woody Wetlands",
        "Emergent Herbaceous Wetlands"
    )
}

