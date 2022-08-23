
library(tidyverse)

d20_raw <- readxl::read_xlsx(
    "data-raw/SC Ticks_Animal Shelters 2020.xlsx", 
    sheet="Color Coded Data",
    na=c("", "N/A")
) |> 
    select(-SC_Region, -EngorgedStatus)

colnames(d20_raw) <- c(
    "date_delivered", "date_received", "shelter_id", "county", 
    "genus", "species", "sex", "life_stage", "count", "comments"
)

d20 <- d20_raw |> 
    drop_na(genus) |> 
    mutate(across(1:2, as.Date), life_stage=str_to_lower(life_stage))

# removing host id because of inconsistent labelling

d21_raw <- readxl::read_xlsx("data-raw/2021_SC Animal Shelter Ticks - USE THIS ONE UPDATED.xlsx") |> 
    select(`Date Picked Up`, `Location / Name`, County, Genus:Number)

colnames(d21_raw) <- c("date_received", "shelter", "county", "genus", "species", "life_stage", "sex", "count")

d21 <- d21_raw |> 
    fill(date_received:county, .direction="down") |> 
    mutate(date_received=as.Date(date_received))

shelters_tidy <- bind_rows(d20, d21) |> 
    mutate(genus=str_to_title(genus), species=str_replace(str_to_lower(species), "unkown", "unknown")) |> 
    relocate(shelter, .after=shelter_id)

write_csv(shelters_tidy, "data-proc/shelters-tidy-20-21.csv")
