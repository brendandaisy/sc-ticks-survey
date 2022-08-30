library(tidyverse)

# Model data preparation----------------------------------------------------------

read_parks_sf <- function(f="geo-files/parks-with-covars.shp", drop=NULL) {
    read_sf(f) |> 
        rename(
            life_stage=lif_stg, land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, 
            min_temp=min_tmp, max_temp=max_tmp, precipitation=prcpttn, jan_min_temp=jn_mn_t
        ) |> 
        # Larvae must absolute be filtered since they were always assigned A. americanum!!!
        filter(species != "maculatum", life_stage != "larva") |> 
        select(-{{drop}})
}

# Initial data preparation--------------------------------------------------------
# in all analyses, A. maculatum, all larvae, and the min_temp---------------------
# covariate are removed-----------------------------------------------------------
prep_parks_model_data <- function(parks_sf) {
    # rescale covariates and format a few things
    parks_rescale <- parks_sf |> 
        mutate(
            across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)), # center and scale
            land_cover=land_cover_labels(land_cover) |> fct_drop(),
            tick_class=ifelse(genus == "Ixodes", str_c(genus, " spp."), str_c(genus, " ", species))
        )
    
    parks_rescale |> 
        group_by(across(c(date, month, site, tick_class, land_cover:mean_rh))) |>
        summarize(pres=ifelse(sum(count) > 0, 1L, 0L), abun=sum(count), .groups="drop") |> 
        mutate(
            id=1:n(), 
            tcnum=case_when(
                tick_class == "Amblyomma americanum" ~ 1L, 
                tick_class == "Dermacentor variabilis" ~ 2L, 
                TRUE ~ 3L
            )
        ) |> 
        arrange(tcnum)
}

# Helpers for working with INLA objects-------------------------------------------

fit_model <- function(formula, data, response="pres", sampling=FALSE, fx_prec=0.3, ...) {
    if (response == "pres") {
        ret <- inla(
            formula,
            data=data,
            family="binomial",
            control.fixed = list(
                expand.factor.strategy="inla",
                prec=fx_prec
            ),
            control.compute=list(config=sampling),
            control.predictor=list(link=1, compute=TRUE), 
            ...
        )
    } else {
        ret <- "not implemented"
    }
    return(ret)
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
        str_replace("Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\.", "intercept")
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
    factor(as.character(lc), levels=lv, labels=lab)
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

