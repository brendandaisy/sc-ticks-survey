# --------------------------------------------------------------------------------
# find-missing-ixodes.R-----------------------------------------------------------
# Calculate the difference between ixodes counts in the parks data,---------------
# and the samples sent to the CDC-------------------------------------------------
# --------------------------------------------------------------------------------

parks <- read_csv("data-proc/parks-data-20-21.csv") # read data as created in tidy-parks-data.R

# Load in the ticks sent to CDC and make it tidy:
prep_cdc_dat <- function(fname) {
    readxl::read_xlsx(fname) |> 
        mutate(date=as.Date(`Date of Collection`), site=`Surveillance Site Name`) |> 
        separate(`Tick Genus species`, c("genus", "species"), sep=" ") |> 
        transmute(
            date, site, genus, species,
            sex=extract_sex(`Specimen Characteristic (for ArboNET)`),
            life_stage=extract_lifestage(`Specimen Characteristic (for ArboNET)`)
        )
}

files <- c(
    "0076_SC21-Ticks_61-111_Field Vial Data.xlsx", 
    "0089_SC21_1-206_Field Vial Data.xlsx", 
    "SC21-Ticks_1-59_Field Vial Data.xlsx"
)

cdc <- map_dfr(files, prep_cdc_dat)

# Now we can actually compare counts
# Matches are made according to date, site, and species
cdc_count <- cdc |> 
    mutate(site=str_remove(site, "( Beach)( State Park)?| State Park")) |> # fix inconsistent naming
    count(date, site, species)

parks_count <- filter(parks, genus == "Ixodes") |> 
    mutate(site=str_remove(site, " Beach| State Park")) |>  # fix inconsistent naming
    count(date, site, species)

diffs <- full_join(cdc_count, parks_count, by=c("date", "site", "species"), suffix=c("_cdc", "_parks")) |> 
    replace_na(list(n_cdc=0, n_parks=0)) |> 
    mutate(diff=n_cdc - n_parks)

write_csv(arrange(diffs, date, site), "ixodes-difference.csv")

ggplot(diffs, aes(interaction(site, date), y=diff, fill=species)) +
    geom_col() +
    scale_x_discrete(guide=guide_axis(angle=45)) +
    labs(y="CDC - Parks", fill="CDC - Parks") +
    theme_bw() +
    theme(plot.margin = margin(2, 2, 2, 90))

ggsave("ixodes-difference.pdf", width=12, height=5.7)
