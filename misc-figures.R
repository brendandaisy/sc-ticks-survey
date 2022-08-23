library(tidyverse)
library(sf)
library(corrplot)
library(car)

parks <- st_read("geo-files/parks-with-covars.shp") |> 
    rename(
        life_stage=lif_stg, land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, 
        min_temp=min_tmp, max_temp=max_tmp, precipitation=prcpttn, jan_min_temp=jn_mn_t
    ) |> 
    filter(species != "maculatum", life_stage != "larva")

parks_rescale <- parks |> 
    mutate(
        across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)), # center and scale
        land_cover=as.factor(land_cover),
        tick_class=ifelse(genus == "Ixodes", str_c(genus, " spp."), str_c(genus, " ", species))
    )

corrplot(cor(select(st_drop_geometry(parks), tree_canopy:mean_rh)), method="number")

parks_model_data <- parks_rescale |> 
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

# VIF with and without land cover
vif(glm(pres ~ land_cover+tree_canopy+elevation+jan_min_temp+min_temp+max_temp+precipitation+mean_rh, "binomial", parks_model_data))
vif(glm(pres ~ tree_canopy+elevation+jan_min_temp+min_temp+max_temp+precipitation+mean_rh, "binomial", parks_model_data))
# Make sure VIF good in both cases with min_temp removed
vif(glm(pres ~ land_cover+tree_canopy+elevation+max_temp+precipitation+mean_rh, "binomial", parks_model_data))
vif(glm(pres ~ tree_canopy+elevation+max_temp+precipitation+mean_rh, "binomial", parks_model_data))

# VIF with and without land cover
vif(glm(abun ~ land_cover+tree_canopy+elevation+jan_min_temp+min_temp+max_temp+precipitation+mean_rh, "poisson", parks_model_data))
vif(glm(abun ~ tree_canopy+elevation+jan_min_temp+min_temp+max_temp+precipitation+mean_rh, "poisson", parks_model_data))
# Make sure VIF good in both cases with min_temp removed
vif(glm(abun ~ land_cover+tree_canopy+elevation+max_temp+precipitation+mean_rh, "poisson", parks_model_data))
vif(glm(abun ~ tree_canopy+elevation+max_temp+precipitation+mean_rh, "poisson", parks_model_data))
