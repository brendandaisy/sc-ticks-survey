library(tidyverse)
library(INLA)
library(corrplot)

occ <- st_read("geo-files/county-model-data.shp")
colnames(occ) <- c(readRDS("geo-files/couny-model-data-colnames.rds"), "geometry")

occ <- occ |> 
    select(-c(life_stage, count_parks, count_shelters)) |> 
    distinct() |> 
    mutate(occ, occ_bin=ifelse(occ == "absent", 0, 1))

pred_mat <- select(st_drop_geometry(occ), land_cover_11:jan_min_temp) |> 
    # select(-c(land_cover_11, land_cover_23, land_cover_41, min_temp, jan_min_temp))
    select(-contains("land_c"), -min_temp, -jan_min_temp)

corrplot(cor(pred_mat), method="number")

terms <- str_c(colnames(pred_mat), collapse="+")

fit1 <- inla(
    as.formula(paste0("occ_bin ~ 0 + tick_class + (", terms, "):tick_class")),
    # occ_bin ~ 0 + tick_class + jan_min_temp:tick_class,
    data=st_drop_geometry(occ),
    family="binomial",
    control.fixed = list(
        expand.factor.strategy="inla",
        prec=0.25
    ),
    control.compute=list(dic=TRUE),
    control.predictor=list(link=1, compute=TRUE)
)

summary(fit1)
occ$occ_bin_pred <- ifelse(fit1$summary.fitted.values$mean > 0.5, 1, 0)

occ_train <- occ |> 
    filter(!is.na(occ_bin))

table(occ_train$occ_bin, occ_train$occ_bin_pred)

fit1$summary.random

ggplot(occ) +
    geom_sf(aes(fill=occ_bin), col="white") +
    geom_sf(aes(fill=occ_bin_pred), data=occ[is.na(occ$occ_bin),], col="orange") +
    facet_wrap(~tick_class) +
    labs(fill="Occupancy status") +
    theme_bw()
