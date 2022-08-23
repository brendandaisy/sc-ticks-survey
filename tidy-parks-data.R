# --------------------------------------------------------------------------------
# tidy-parks-data.R---------------------------------------------------------------
# Initial cleaning and combining of the '20 and '21 data into a consistent format-
# --------------------------------------------------------------------------------

# Load required packages
library(tidyverse)

# Helper functions----------------------------------------------------------------

# Helpers for 2021
extract_sex <- function(x) {
    ifelse(
        str_detect(x, "[Ff]emale"), "F",
        ifelse(str_detect(x, "[Mm]ale"), "M", "I")
    )
}

extract_lifestage <- function(x) {
    ifelse(
        str_detect(x, "[Aa]dult"), "adult", 
        ifelse(str_detect(x, "[Nn]ymph"), "nymph", "larva")
    )
}

extract_samp_method <- function(x, y) {
    map2_chr(x, y, ~{
        if (!is.na(.y))
            if (str_detect(.y, "[bB]ody")) "body" else "bite"
        else
            if (str_detect(.x, "[dD]rag")) "drag" else "CO2"
    })
}

d20_locations <- tribble(
    ~site, ~long, ~lat,
    "Wannamaker", 32.97732072931062, -80.0542338023381,
    "Hunting Island", 32.366218592365534, -80.44883968623985,
    "McAlhany Nature Preserve", 33.15216929240691, -80.6952231079484,
    "Paris Mountain", 34.941541558490535, -82.38950431363385,
    "Woods Bay", 33.95079112665409, -79.97337555978886,
    "Goff Swamp Hunting Club", 34.02762373255747, -80.68091375523846
)

d21_raw <- readxl::read_xlsx("data-raw/2021 SC Ticks Parks CDC Format.xlsx", col_names=d21_cols, skip=2) |> 
    drop_na(site)

d21_locations <- d21_raw |> 
    drop_na(lat, long) |> 
    distinct(site, lat, long) |> # to give each entry equal weight when calc. mean
    group_by(site) |>
    summarise(lat=mean(lat), long=mean(long))

locations <- bind_rows(d20_locations, d21_locations)

# Clean the 2020 data-------------------------------------------------------------
d20_cols <- c(
    "date", "epi_week", "collector", "site", "county", "data_entry", "confirmed", "temp_min",
    "temp_max", "rh_pct", "precipitation", "elevation"
)

d20_tick_cols <- readxl::read_xlsx("data-raw/SC Ticks 2020.xlsx", sheet="TOTAL Tick Collections") |> 
    select(Drag_Genus:Bite_Number) |> 
    colnames()

d20_cols <- c(d20_cols, d20_tick_cols, "comments")

d20 <- readxl::read_xlsx("data-raw/SC Ticks 2020.xlsx", sheet="TOTAL Tick Collections", col_names=d20_cols, skip=1) |> 
    drop_na(date) |>  # drop the empty rows at the bottom
    mutate(
        date=as.Date(date), # convert to better date fmt
        temp_min=as.character(temp_min),
        rh_pct=as.character(rh_pct), precipitation=as.character(precipitation)
    )

drag <- d20 |> 
    select(date:elevation, contains("Drag"), comments) |> 
    drop_na(Drag_Genus) |> # this assumes all >0 counts have been identified
    rename(genus=Drag_Genus, species=Drag_Species, sex=Drag_Sex, life_stage=Drag_LifeStage, count=Drag_Number) |> 
    mutate(sampling_method="drag")

co2 <- d20 |> 
    select(date:elevation, contains("CO2"), comments) |> 
    # drop_na(CO2_Genus) |> 
    rename(genus=CO2_Genus, species=CO2_Species, sex=CO2_Sex, life_stage=CO2_LifeStage, count=CO2_Number) |> 
    mutate(sampling_method="CO2")

body <- d20 |> 
    select(date:elevation, contains("Body"), comments) |> 
    drop_na(Body_Genus) |> 
    rename(genus=Body_Genus, species=Body_Species, sex=Body_Sex, life_stage=Body_LifeStage, count=Body_Number) |> 
    mutate(sampling_method="body")

bite <- d20 |> 
    select(date:elevation, contains("Bite"), comments) |> 
    drop_na(Bite_Genus) |> 
    rename(genus=Bite_Genus, species=Bite_Species, sex=Bite_Sex, life_stage=Bite_LifeStage, count=Bite_Number) |> 
    mutate(sampling_method="bite")

d20_tidy <- bind_rows(drag, co2, body, bite) |> 
    mutate(site=ifelse(site == "Edisto Beach", "Edisto", site)) |> 
    left_join(locations)

# Do the 2021 cleaning------------------------------------------------------------
d21_cols <- c(
    "site", "street", "city", "state", "county", "epi_week",
    "collector", "time_of_day", "temp_min", "temp_max", "rh_pct",
    "precipitation", "data_entry", "confirmed", "lat", "long", "habitat",
    "ownership", "method", "method_extra", "date", "genus", "species", "adult_female",
    "adult_male", "adult_unknown", "nymph", "larvae_collected", "larvae_counted", "larvae",
    "sampling_notes", "dist_covered", "cdc_id", "comments"
)

d21_tidy <- d21_raw |>
    mutate(date=as.Date(date)) |> 
    pivot_longer(c(adult_female:nymph, larvae), values_to="count") |> 
    # filter(!is.na(genus), name != "larvae" | (name == "larvae" & larvae_counted == "Yes")) |> 
    # filter(!is.na(genus)) |> 
    mutate(
        sex=extract_sex(name), 
        life_stage=extract_lifestage(name), 
        sampling_method=extract_samp_method(method, method_extra)
    ) |> 
    select(!c(name, method, method_extra, lat, long)) |> 
    left_join(locations) # add the previously computed long lat (from means of orig entries)

# Save both data sets into one dataframe------------------------------------------

parks <- bind_rows(d20_tidy, d21_tidy) |> 
    select(date, site, intersect(colnames(d20_tidy), colnames(d21_tidy)), everything())

write_csv(parks, "data-proc/parks-data-20-21.csv")
