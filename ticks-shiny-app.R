library(tidyverse)
library(sf)
library(shiny)
library(plotly)

source("shiny-app-helpers.R")

# Global Variables---------------------------------------------------------
parks <- read_parks_data("data-raw/2021 SC Ticks Parks CDC Format.xlsx")

sc_counties <- st_read("geo-files/south-carolina-county-boundaries.shp") |>
    st_transform("+proj=longlat +datum=WGS84")

park_locations <- parks |>
    select(site, lat, long) |>
    distinct(site, .keep_all=TRUE) |>
    st_as_sf(crs="+proj=longlat +datum=WGS84", coords=c("lat", "long")) |>
    st_transform(st_crs(sc_counties))

# User Interface-----------------------------------------------------------
js_script <- '
    var dimension = [0, 0];
    $(document).on("shiny:connected", function(e) {
    dimension[0] = window.innerWidth;
    dimension[1] = window.innerHeight;
    Shiny.onInputChange("dimension", dimension);
    });
    $(window).resize(function(e) {
    dimension[0] = window.innerWidth;
    dimension[1] = window.innerHeight;
    Shiny.onInputChange("dimension", dimension);
    });
'

ui <- fluidPage(
    title = "2021 South Carolina Tick Survey",
    tags$head(tags$script(js_script)),
    sidebarLayout(
        mainPanel(plotOutput("plot_obj", width="100%"), width=9), # may have to change height dep on screen
        sidebarPanel(
            selectInput("tick_col" , h3("Tick Column"), colnames(select(parks, adult_female:nymph))),
            downloadButton("downloadPlot", "Download Plot"),
            width=3
        )
    )
)

get_months_totals <- function(tick_col) {
    parks |>
        group_by(month=fct_inorder(format(date, "%B")), site=fct_inorder(site), .drop=FALSE) |>
        summarise(lat, long, total := sum( .data[[tick_col]] )) |>
        ungroup() |>
        distinct(month, site, .keep_all=TRUE) |>
        st_as_sf(crs="+proj=longlat +datum=WGS84", coords=c("lat", "long")) |>
        st_transform(st_crs(sc_counties))
}


# month_totals <- parks |>
#     group_by(month=fct_inorder(format(date, "%B")), site=fct_inorder(site), .drop=FALSE) |>
#     summarise(lat, long, adult_total=sum(adult_female + adult_male + adult_unknown)) |>
#     ungroup() |>
#     distinct(month, site, .keep_all=TRUE)

# exp <- expand(tmp, month, nesting(site, lat, long)) |>
#     distinct(month, site, .keep_all=TRUE)

# month_totals <- anti_join(exp, month_totals) |>
#     mutate(adult_total=NA) |>
#     bind_rows(tmp) |>
#     st_as_sf(crs="+proj=longlat +datum=WGS84", coords=c("lat", "long")) |>
#     st_transform(st_crs(sc_counties))

match_county_month <- function(m, ct_idx, mems, mtot) {
    matches <- filter(mtot, mems[,ct_idx], month == m)$total
    if (length(matches) == 0)
        return(NA)
    return(sum(matches))
}

plot_county_cts_month <- function(tick_col) {
    month_totals <- get_months_totals(tick_col)

    county_month_long <- sc_counties |>
        mutate(month=list(unique(month_totals$month))) |>
        unnest(month)

    mems <- st_within(month_totals, county_month_long, sparse=FALSE)
    match_func <- function(m, idx) match_county_month(m, idx, mems, month_totals)

    county_month_long <- county_month_long |>
        mutate({{tick_col}} := imap_dbl(month, match_func))

    gg <- ggplot() +
        geom_sf(aes(fill=log1p(.data[[tick_col]])), data=county_month_long) +
        theme_bw() +
        facet_wrap(~month, nrow=3) +
        labs(fill=paste0(str_to_title(str_replace(tick_col, "_", " ")), " (log scale)"))
    
    # geom_sf_text(aes(label=site), data=parks_proj, size=4, col="gray90") +
    return(gg)
}

# plot_county_cts_month("adult_female")

# Do work for current input and return plots
server <- function(input, output) {
    observeEvent(input$dimension,{
        output$plot_obj <- renderPlot({ # gets called when inputs are changed
            gg <- plot_county_cts_month(input$tick_col)
            gg
        }, width=0.7*as.numeric(input$dimension[1]), height=0.7*as.numeric(input$dimension[2])
    )})

    output$downloadPlot <- downloadHandler(
        filename = function() paste0("county_", input$tick_col, ".pdf"), 
        content = function(file) ggsave(file, plot_county_cts_month(input$tick_col))
    )
}

# Run the app
shinyApp(ui = ui, server = server)
