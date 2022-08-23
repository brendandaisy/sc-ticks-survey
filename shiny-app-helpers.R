
read_parks_data <- function(fname, fill_geo=TRUE) {
    cols <- c(
        "site", "street", "city", "state", "county", "epi_week",
        "collector", "time_of_day", "temp_min", "temp_max", "rh_pct",
        "precipitation", "data_entry", "confirmed", "lat", "long", "habitat",
        "ownership", "method", "method_extra", "date", "genus", "species", "adult_female",
        "adult_male", "adult_unknown", "nymph", "larvae_collected", "larvae_counted", "larvae",
        "sampling_notes", "dist_covered", "cdc_id", "comments"
    )

    ret <- readxl::read_xlsx(fname, col_names=cols, skip=2) |>
        drop_na(site) |>
        mutate(date=as.Date(date))
    
    if (fill_geo)
        ret <- fill_missing_coords(ret)
    return(ret)
}

fill_missing_coords <- function(df) {
    uniques <- df |>
        select(site, lat, long) |>
        distinct() |>
        drop_na(lat, long)
    
    # replace geocoords that are NA with the unique site associated with those coords:
    na_geo <- which(is.na(df$lat) | is.na(df$long))
    na_site <- df$site[na_geo]
    df[na_geo, c("lat", "long")] <- uniques[match(na_site, uniques$site), c("lat", "long")]
    
    return(df)
}
