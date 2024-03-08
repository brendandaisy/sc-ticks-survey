library(tidyverse)
library(sf)

# Summary/overview of initial collections (not gonna use)

read_csv("data-proc/parks-data-20-21.csv") |> 
    filter(life_stage != "larva") |> # must be filtered since all were assigned A americanum
    mutate(
        month=lubridate::month(month, label=TRUE, abbr=FALSE),
        year=as.factor(year),
        species=str_replace(species, "affinis", "kearnsi"), # new taxonomy
        tick_class=str_c(genus, " ", species),
        pres=ifelse(count > 0, 1, 0)
    ) |> 
    select(-c(genus, species, count))
    

parks_obs <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) |> 
    mutate(
        month=lubridate::month(month, label=TRUE, abbr=FALSE),
        year=as.factor(lubridate::year(date))
        # pres=fct_recode(as.character(pres), "absent"="0", "present"="1")
    )

visits <- distinct(parks_obs, date, site, month, year)

visits |> 
    count(month, year) |> 
    ggplot(aes(month, n, fill=year)) +
    geom_col(position=position_stack()) +
    theme_bw()

sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(st_crs(parks_obs)) |> 
    st_union() |> 
    nngeo::st_remove_holes()

parks_obs |> 
    count(site) |> 
    ggplot(aes(size=n, col=n)) +
    geom_sf(data=sc_state, fill=NA, col="gray70", linewidth=0.7, alpha=0.9) +
    geom_sf() +
    scale_color_viridis_c(, option="mako") +
    labs(col="count", title="total number of visits") +
    scale_size(guide="none") +
    theme_bw() +
    theme(
        # legend.position="top",
        panel.grid=element_blank(),
        axis.text=element_blank(), axis.ticks=element_blank(),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
    )

psumm <- parks_obs |> 
    group_by(month, site, tick_class) |> 
    summarise(pres=ifelse(any(pres == 1), "present", "absent"), .groups="drop")

sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(st_crs(parks_obs)) |> 
    st_union() |> 
    nngeo::st_remove_holes()

ggplot(psumm, aes(col=pres)) +
    geom_sf(data=sc_state, fill=NA, col="gray70", linewidth=0.7, alpha=0.9) +
    geom_sf(size=1.1) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    scale_color_manual(values=c("#3798a9", "#f53dbc")) +
    labs(col="Presence of nymphs during visit") +
    theme_bw() +
    theme(
        strip.text=element_text(size=rel(1.05)),
        axis.text=element_blank(), axis.ticks=element_blank(),
        panel.spacing=unit(0, "mm"),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
        legend.position="top",
        legend.box.margin=unit(c(0, 0, 0, 0), "mm")
    )

### OLD STUFF
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

# For ASTMH Poster----------------------------------------------------------------

df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

all_loc_sp <- read_sf("data-proc/parks-design-space.shp") |> 
    distinct(site, geometry)

risk_mod_sp <- semi_join(all_data, risk_mod, by=c("month", "site")) |> 
    mutate(month=lubridate::month(date, label=TRUE, abbr=TRUE))

sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(st_crs(all_loc_sp)) |> 
    st_union() |> 
    nngeo::st_remove_holes()

library(ggthemes)

ggplot(risk_mod_sp) +
    geom_sf(data=sc_state, fill=NA, col="#448dc6ff", linewidth=0.9, alpha=0.9) +
    geom_sf(col="white", size=2) +
    # geom_sf(col="#ff856dff", size=2.5, shape=1) +
    geom_sf(col="#ff856dff", size=1) +
    theme_map() +
    theme(axis.text=element_blank(), axis.ticks=element_blank())

ggsave("astmh-poster/risk-risk.pdf", width=2, height=2)
