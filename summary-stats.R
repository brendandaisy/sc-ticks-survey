# summary-stats.R-----------------------------------------------------------------
# Simple summary statistics about the 2020-21 parks data, reported in the text----
# --------------------------------------------------------------------------------

source("other-helpers.R")

parks_tidy <- read_csv("data-proc/parks-data-20-21.csv")
parks_obs <- read_parks_sf("data-proc/parks-observed.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE)


vis <- distinct(parks_tidy, date, site, .keep_all=TRUE)  # keeps the first row of values for each visit
    
# number of visits which happened each year
count(vis, year)

# number/percentage of visits each tick was found
parks_tidy |> 
    filter(count > 0) |> 
    count(genus, species)

parks_obs |> 
    filter(pres > 0) |> 
    count(tick_class) |> 
    mutate(p=100*n/nrow(vis))

    