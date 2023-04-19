library(tidyverse)

# Initial data preparation--------------------------------------------------------
# in all analyses, A. maculatum, all larvae, and the min_temp---------------------
# covariate are removed-----------------------------------------------------------
read_parks_sf <- function(f="data-proc/parks-observed.shp", drop=NULL) {
    read_sf(f) |> 
        rename(
            life_stage=lif_stg, land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, 
            min_temp=min_tmp, max_temp=max_tmp, precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        # Larvae must absolute be filtered since they were always assigned A. americanum!!!
        filter(species != "maculatum", life_stage != "larva") |> 
        relocate(land_cover, tree_canopy, elevation, min_temp, max_temp, precipitation, jan_min_temp, mean_rh, .after=life_stage) |> 
        select(-{{drop}})
}

prep_parks_model_data <- function(parks_sf, rescale=TRUE) {
    # rescale covariates and format a few things
    parks_rescale <- parks_sf |> 
        mutate(
            across(tree_canopy:mean_rh, ~if (rescale) (.x - mean(.x)) / sd(.x) else .x), # center and scale
            land_cover=land_cover_labels(land_cover),
            tick_class=ifelse(genus == "Ixodes", str_c(genus, " spp."), str_c(genus, " ", species))
        )
    
    ret <- parks_rescale |> 
        group_by(across(c(date, month, site, tick_class, land_cover:mean_rh))) |>
        summarize(pres=ifelse(sum(count) > 0, 1L, 0L), abun=sum(count), .groups="drop")
    
    mutate(ret, tcnum=case_when(
        tick_class == "Amblyomma americanum" ~ 1L, 
        tick_class == "Dermacentor variabilis" ~ 2L, 
        TRUE ~ 3L
    )) |> 
        arrange(tcnum)
}

append_pred_data <- function(obs_df, f="data-proc/parks-design-space.shp", pred_grid=NULL, drop_new_lcc=TRUE) {
    pred_df <- st_read(f) |> 
        rename(
            land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
            precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        mutate(land_cover=land_cover_labels(land_cover)) |> 
        relocate(land_cover, tree_canopy, elevation, min_temp, max_temp, precipitation, jan_min_temp, mean_rh, .after=date) |> 
        select(-min_temp)
    
    # Add a unique site label to each grid location (i.e. const. over time)
    pred_df <- pred_df |> 
        mutate(data=list(tibble(tick_class=unique(obs_df$tick_class), pres=NA, tcnum=1:3))) |> 
        unnest(c(data))
    
    # Optionally add the pred grid. WARNING: this will result in different scaling than if pred grid were not included
    if (!is.null(pred_grid))
        pred_df <- bind_rows(pred_df, pred_grid)
    
    if (drop_new_lcc)
        pred_df <- filter(pred_df, land_cover %in% levels(obs_df$land_cover)) |> mutate(land_cover=fct_drop(land_cover))
    
    # WARNING: all_data has now been scaled, and is the reference covar values for all subdesigns
    return(prep_new_data(obs_df, pred_df, scale=TRUE))
}

# parks_data should already be centered/scaled (which is fine because this data is never actually influencing model fit),
# and should contain both the visitation and design space parks data
prep_pred_grid <- function(parks_data, f="geo-files/covar-grid-16km.shp") {
    covar_grid <- st_read(f) |> 
        rename(
            land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, min_temp=min_tmp, max_temp=max_tmp,
            precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        mutate(land_cover=land_cover_labels(land_cover)) |> 
        select(-min_temp)
    
    # Add a unique site label to each grid location (i.e. const. over time)
    covar_grid <- covar_grid |> 
        group_by(geo=as.character(geometry)) |> 
        mutate(site=str_c("g", cur_group_id())) |> 
        ungroup() |> 
        select(-geo) |> 
        mutate(data=list(tibble(tick_class=unique(parks_data$tick_class), pres=NA, tcnum=1:3))) |> 
        unnest(c(data))
    
    covar_grid
}

# This should be a function, since 1) the covariate scaling is very much dependent on the chosen locations now, and 
# 2) the joint residuals require the data to be arranged each time
prep_new_data <- function(known_df, new_df=NULL, scale=FALSE) {
    ret <- bind_rows(known_df, new_df)
    
    if (scale)
        ret <- mutate(ret, across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)))
    return(ret)
}

# Helpers for working with INLA objects-------------------------------------------

fit_model <- function(formula, data, fx_prec, response="pres", control_compute=list(mlik=FALSE), ...) {
    if (is(data, "sf"))
        data <- st_drop_geometry(data)
    if (response == "pres") {
        ret <- inla(
            formula,
            data=data,
            family="binomial",
            control.fixed = list(
                expand.factor.strategy="inla", # nec. since alternative removes unused lvls
                prec=fx_prec
            ),
            control.compute=control_compute,
            control.predictor=list(link=1, compute=TRUE), 
            ...
        )
    } else {
        ret <- "not implemented"
    }
    return(ret)
}

# optional update formula for random effects
formula_jsdm <- function(df, rx=~.) {
    ret <- pres ~ -1 + tick_class +  
        tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation +
        tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh +
        f(id, model="iid3d", n=nrow(df), hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3)))))
    
    update(ret, rx)
}

# the (current...) preferred model for performing BED
formula_bed <- function(...) {
    prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.1)))

    pres ~ -1 + tick_class + land_cover + tree_canopy + elevation + max_temp + 
        precipitation + jan_min_temp + mean_rh + 
        f(month, model = "ar1", hyper = prec_pri, group = tcnum, control.group = list(model = "iid", hyper = prec_pri)) + 
        f(site, model = "iid", hyper = prec_pri, group = tcnum, control.group = list(model = "iid", hyper = prec_pri))
}

fit_all_data <- function(all_data) {
    f <- formula_bed()
    fit_model(f, all_data, fx_prec=0.2, control_compute=list(mlik=TRUE, dic=TRUE))
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
    lab <- c(
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
    ret <- factor(as.character(lc), levels=lv, labels=lab)
    return(fct_drop(ret))
}

land_cover_palette <- function() {
    c(
        "Open Water" = "#5475a8",
        "Perennial Ice/Snow" = "white",
        "Developed, Open Space" = "#e8d1d1",
        "Developed, Low Intensity" = "#e29e8c",
        "Developed, Medium Intensity" = "#ff0000",
        "Developed, High Intensity" = "#b50000",
        "Barren Land" = "#d2cdc0",
        "Deciduous Forest" = "#85c77e",
        "Evergreen Forest" = "#38814e",
        "Mixed Forest" = "#d4e7b0",
        "Dwarf Scrub" = "#af963c",
        "Shrub/Scrub" = "#dcca8f",
        "Grassland/Herbaceous" = "#fde9aa",
        "Sedge/Herbaceous" = "#d1d182",
        "Lichens" = "#a3cc51",
        "Moss" = "#82ba9e",
        "Pasture/Hay" = "#fbf65d",
        "Cultivated Crops" = "#ca9146",
        "Woody Wetlands" = "#c8e6f8",
        "Emergent Herbaceous Wetlands" = "#64b3d5"
    )
}

