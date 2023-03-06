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

prep_parks_model_data <- function(parks_sf, rescale=TRUE, joint_compat=TRUE) {
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
    
    if (joint_compat) {
        ret <- mutate(ret, tcnum=case_when(
            tick_class == "Amblyomma americanum" ~ 1L, 
            tick_class == "Dermacentor variabilis" ~ 2L, 
            TRUE ~ 3L
        )) |> 
            arrange(tcnum) |> 
            mutate(id=1:n())
    }
    return(ret)
}

append_pred_grid <- function(parks_df, f="geo-files/covar-grid.shp") {
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
        mutate(data=list(tibble(tick_class=unique(parks_df$tick_class), pres=NA, tcnum=1:3))) |> 
        unnest(c(data))
    
    # Add levels missing from parks data but present in SC
    levels(parks_df$land_cover) <- levels(covar_grid$land_cover)
    
    # Predict expected risk across grid, for sampling---------------------------------
    
    # WARNING: all_data has now been scaled, and is the reference covar values for all subdesigns
    return(prep_new_data(parks_df, covar_grid, scale=TRUE))
}

# This should be a function, since 1) the covariate scaling is very much dependent on the chosen locations now, and 
# 2) the joint residuals require the data to be arranged each time
prep_new_data <- function(known_df, new_df=NULL, scale=FALSE) {
    ret <- known_df |>
        bind_rows(new_df) |>
        arrange(tcnum) |>
        mutate(id=1:n())
    
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
    pres ~ -1 + tick_class +  
        tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation +
        tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh
}

get_fixed_effects <- function(..., new_loc) {
    amb <- get("tick_classAmblyomma americanum") +
        get(paste0("tick_classAmblyomma americanum:land_cover", new_loc$land_cover)) +
        get("tick_classAmblyomma americanum:tree_canopy")*new_loc$tree_canopy +
        get("tick_classAmblyomma americanum:jan_min_temp")*new_loc$jan_min_temp +
        get("tick_classAmblyomma americanum:max_temp")*new_loc$max_temp +
        get("tick_classAmblyomma americanum:precipitation")*new_loc$precipitation +
        get("tick_classAmblyomma americanum:mean_rh")*new_loc$mean_rh
    
    der <- get("tick_classDermacentor variabilis") +
        get(paste0("tick_classDermacentor variabilis:land_cover", new_loc$land_cover)) +
        get("tick_classDermacentor variabilis:tree_canopy")*new_loc$tree_canopy +
        get("tick_classDermacentor variabilis:jan_min_temp")*new_loc$jan_min_temp +
        get("tick_classDermacentor variabilis:max_temp")*new_loc$max_temp +
        get("tick_classDermacentor variabilis:precipitation")*new_loc$precipitation +
        get("tick_classDermacentor variabilis:mean_rh")*new_loc$mean_rh
    
    ix <- get("tick_classIxodes spp.") +
        get(paste0("tick_classIxodes spp.:land_cover", new_loc$land_cover)) +
        get("tick_classIxodes spp.:tree_canopy")*new_loc$tree_canopy +
        get("tick_classIxodes spp.:jan_min_temp")*new_loc$jan_min_temp +
        get("tick_classIxodes spp.:max_temp")*new_loc$max_temp +
        get("tick_classIxodes spp.:precipitation")*new_loc$precipitation +
        get("tick_classIxodes spp.:mean_rh")*new_loc$mean_rh
    
    return(c(amb, der, ix))
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
        str_remove("tick_class") |> 
        str_remove("(Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\.):") |> 
        str_replace("Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\.", "intercept") |> 
        str_remove("land_cover") |> 
        fct_relevel("intercept", after=Inf) |> 
        fct_relevel("tree_canopy", "elevation", "jan_min_temp", "max_temp", "precipitation", "mean_rh", after=0)
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

