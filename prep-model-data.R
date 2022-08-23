# --------------------------------------------------------------------------------
# prep-model-data.R---------------------------------------------------------------
# Take tidy parks data and add covariates, drop some not useful columns, remove---
# Dermacentor observations, and convert implicit zero-counts to explicit----------
# --------------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(raster, exclude=c("select"))
library(furrr)
library(lubridate)

# Helper functions----------------------------------------------------------------
# project_parks <- function(parks_tidy, ref_raster="mean-temp", pattern="bil$") {
#     ret <- parks_tidy |>
#         select(site, lat, long) |>
#         distinct(site, .keep_all=TRUE) |>
#         st_as_sf(crs="+proj=longlat +datum=WGS84", coords=c("lat", "long"))
#     
#     files <- list.files(paste0("geo-files/", ref_raster), pattern=pattern, full.names=TRUE)
#     return(st_transform(ret, crs(raster(files[1]))))
# }

agg_parks <- function(parks_tidy, groups_inc=c("Ixodes_adult", "Amblyomma_adult", "Amblyomma_nymph")) {
    parks_tidy |> 
        select(date, site, genus:count, sampling_method, lat, long) |> 
        mutate(month=month(floor_date(date, "months"))) |> 
        unite("group", genus, life_stage) |> 
        filter(group %in% groups_inc) |> # A. larva too since very few collection dates
        group_by(date, site, group, month, lat, long) |> # marginalize over factors not considered 
        summarise(count=sum(count)) |> 
        ungroup() |> 
        complete(nesting(date, site, month, lat, long), group, fill=list(count=0))
}

project_parks <- function(parks, ref_raster="mean-temp", pattern="bil$") {
    files <- list.files(paste0("geo-files/", ref_raster), pattern=pattern, full.names=TRUE)
    
    parks |>
        st_as_sf(crs="+proj=longlat +datum=WGS84", coords=c("lat", "long")) |>
        st_transform(crs(raster(files[1])))
}

project_shp <- function(shp_file, ref_raster="mean-temp", pattern="bil$") {
    ret <- st_read(shp_file)
    files <- list.files(paste0("geo-files/", ref_raster), pattern=pattern, full.names=TRUE)
    return(st_transform(ret, crs(raster(files[1]))))
}

get_landcover_grid <- function(bounds, saved_file=NULL) {
    if (!is.null(saved_file))
        return(st_read(saved_file))
    
    lcov <- raster("geo-files/land-cover/NLCD_2019_Land_Cover_L48_20210604_yDj0p7IM7PAX2WUUipna.tiff")
    lcov_4k <- aggregate(lcov, 133, fun=modal)
    
    lcov_sp <- rasterToPoints(lcov_4k, spatial=TRUE) |> 
        st_as_sf() |> 
        rename(land_cover=NLCD_2019_Land_Cover_L48_20210604_yDj0p7IM7PAX2WUUipna) |> 
        st_transform(crs(bounds))
    
    u <- st_within(lcov_sp, bounds)
    return(lcov_sp[map_lgl(u, ~length(.x) > 0),])
}

get_raster_stack <- function(var, dates) {
    file_dates <- str_replace_all(dates, "-", "")
    pattern <- str_c(".+_(", str_c(file_dates, collapse="|"), ").+\\.bil$")
    files <- list.files(paste0("geo-files/", var), pattern=pattern, full.names=TRUE)
    as.list(stack(files))
}

load_prism_var <- function(var, dates, nlcd_point_grid) {
    stack <- get_raster_stack(var, dates)
    varname <- str_replace_all(var, "-", "_")
    
    furrr::future_map2_dfr(
        stack, dates, 
        extract_raster_grid, nlcd_point_grid, 
        .options=furrr_options(seed=NULL)
    ) |> 
        rename({{varname}} := rval)
}

# load_prism_var_park <- function(park_sp, var) {
#     stack <- get_raster_stack(var, unique(park_sp$date))
#     varname <- str_replace_all(var, "-", "_")
#     
#     park_sp <- group_by(park_sp, date) |> 
#         mutate(date_ind=cur_group_id()) |> 
#         ungroup()
# 
#     rvals <- map2_dbl(
#         park_sp$date_ind, park_sp$geometry, ~extract_raster_park(stack[[.x]], .y),
#         .options=furrr_options(seed=NULL)
#     )
# 
#     mutate(park_sp, {{varname}} := rvals)
# }
# 
# extract_raster_park <- function(r, point) {
#     raster::extract(r, as(point, "Spatial"))
# }

extract_raster_grid <- function(r, date, point_grid) {
    rvals <- raster::extract(r, point_grid)
    
    point_grid |> 
        st_as_sf() |> 
        mutate(rval=rvals, date=date)
}

relative_humidity <- function(temp, dew, beta=17.625, lambda=243.04) {
    num <- exp(beta * dew / (lambda + dew))
    denom <- exp(beta * temp / (lambda + temp))
    100 * num / denom
}

## Plan: prep "observed" and "unobserved" points seperately
## unobserved points consist of grid over state, with meas. taken once at begin. of month
## Plan: prep the grid with month/group/count but not date
## Send with dates of each month to extract full grid for each data
## appending rows then give all info in grid so we can just
## bind with the site/date observations

plan(multisession, workers=4) # set up parallel session

parks <- read_csv("data-proc/parks-data-20-21.csv") # read data as created in tidy-parks-data.R
parks_agg <- agg_parks(parks) |> 
    filter(month(date) <= 10)
# park_locations <- project_parks(parks, "mean-temp") # now can get reference CRS from here
parks_proj <- project_parks(parks_agg, "mean-temp") # now can get reference CRS from here

# Update parks points with the matching land_cover value (easier now than later)
parks_proj <- st_join(parks_proj, lcov_grid, st_nearest_feature)

dates_parks <- distinct(parks_proj, date) |> slice(1:(n()-2)) |> pull()
park_locations <- distinct(parks_proj, site, geometry)
mtemp_park <- load_prism_var("mean-temp", dates_parks, park_locations) # No pipe cause each line expensive!
precip_park <- load_prism_var("precipitation", dates_parks, park_locations)
dew_park <- load_prism_var("dew-point", dates_parks, park_locations)

covars_park <- mtemp_park |> 
    mutate(precipitation=precip_park$precipitation, mean_rh=relative_humidity(mean_temp, dew_park$dew_point))

covars_park
    
sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(crs(parks_proj)) |> 
    st_union() |> 
    nngeo::st_remove_holes() # remove a few holes where counties didn't touch

## If grid not saved:
# lcov_grid <- get_landcover_grid(sc_state) # grid will be proj to state (which should be proj to parks_proj)
# st_write(lcov_grid, "geo-files/landcover-grid.shp")
## Otherwise:
lcov_grid <- get_landcover_grid(NULL, saved_file="geo-files/landcover-grid.shp")


# Update grid points with PRISM data
dates_grid <- parks |> 
    pull(date) |> 
    floor_date("months") |> 
    unique()

lcov_grid <- lcov_grid |> 
    mutate(site=NA, group=list(c("Ixodes_adult", "Amblyomma_adult", "Amblyomma_nymph")), count=NA)

# covar_pts <- bind_rows(lc_grid_proj, st_join(park_locations, lc_grid_proj, st_nearest_feature))
grid_sp <- as(lcov_grid, "Spatial")
mtemp <- load_prism_var("mean-temp", dates_grid[1:(length(dates_grid)-2)], grid_sp) # No pipe cause each line expensive!
precip <- load_prism_var("precipitation", dates_grid[1:(length(dates_grid)-2)], grid_sp)
dew <- load_prism_var("dew-point", dates_grid[1:(length(dates_grid)-2)], grid_sp)

ggplot(filter(precip, date == "2020-04-01")) +
    geom_sf(aes(col=precipitation))

covar_grid <- mtemp |>
    mutate(precipitation=precip$precipitation, mean_rh=relative_humidity(mean_temp, dew$dew_point)) |> 
    unnest(c(group))

## NOW, extract exclusively the observed date/site pairs

st_write(covar_grid, "geo-files/TODO.shp")

## Finally make the full dataset!!
full_model_data <- parks_proj |> 
    left_join(st_drop_geometry(covars_park)) |> 
    bind_rows(covar_grid)

st_write(full_model_data, "geo-files/full-model-data.shp")

# get_raster_var <- function(var, locations, pattern="bil$") {
#     files <- list.files(paste0("geo-files/", var), pattern=pattern, full.names=TRUE)
#     varstack <- as.list(stack(files))
#     varname <- str_replace(var, "-", "_")
#     locations <- as(st_transform(locations, crs(varstack[[1]])), "Spatial") # make sure CRS match and revert to sp
#     
#     furrr:::future_map_dfr(varstack, rextract_points, locations, .options=furrr_options(seed=NULL)) |> 
#         rename({{varname}} := vals)
# }

# rextract_points <- function(r, locations) {
#     vals <- raster::extract(r, locations)
#     date <- as.Date(str_extract(r@file@name, "\\d{8}"), "%Y%m%d")
#     tibble(site=locations$site, vals=vals, date=date)
# }



## Add environmental covariates
plan(multisession, workers=4) # set up parallel session

mtemp <- get_raster_var("mean-temp", park_locations)
precip <- get_raster_var("precipitation", park_locations)
dew <- get_raster_var("dew-point", park_locations)

covars <- reduce(list(mtemp, precip, dew), left_join) |> 
    mutate(mean_rh=relative_humidity(mean_temp, dew_point))



lcov_sp <- rasterToPoints(lcov_4k) |> 
    st_as_sf() |>
    rename(land_cover=NLCD_2019_Land_Cover_L48_20210604_yDj0p7IM7PAX2WUUipna)

sc_counties <- 



ggplot(lcov_sp) +
    geom_sf(aes(fill=as.factor(land_cover), col=as.factor(land_cover))) +
    geom_sf(data=sc_state2, fill=NA)

## Match daily covars with park data
parks_mod <- parks |> 
    select(date, site, genus:count, sampling_method) |> 
    mutate(month=month(floor_date(date, "months"))) |> 
    filter(genus != "Dermacentor") |> 
    complete(nesting(date, site, species, sex, sampling_method, month), genus, life_stage, fill=list(count=0)) |> 
    left_join(covars, by=c("site", "date")) |> # add in daily covariates
    drop_na(mean_temp:mean_rh) |>   # this will remove the nov and dec collections since no available covars ATM
    select(date, month, site, genus, species, life_stage, count, everything())
    
write_csv(parks_mod, "data-proc/parks-model-data.csv")

## Pick up where we left off...add the ecoregions
parks_mod <- read_csv("data-proc/parks-model-data.csv")
ecoreg <- project_shp("geo-files/ecoregion/sc_eco_l4.shp", "mean-temp")

parks_eco <- park_locations |> 
    mutate(ecoreg = map_chr(geometry, ~ecoreg$US_L4NAME[st_within(.x, ecoreg, sparse=FALSE)])) |> 
    st_drop_geometry() # don't need to merge it in to the data

write_csv(left_join(parks_mod, parks_eco), "data-proc/parks-model-data.csv")

## tidy to CDC script

## email people about tagging along in lab

## maddie and kia (talk to maddie)

## add veternary data
# not close to parks
# spatial points seems like way to go
# county level vet data
# double check: one line per one dog? counts
# only check when dogs come in
# = peridomestic/domestic environment
# human population density (block/block groups/zip code); county at first?
# ^census data: geolytics (2010 boundaries)
# engorgement = pseudotime. eng=3-10days, part=2, un=1 (Google it!)
# on later dogs, no. ticks were collected per dog, but before []

