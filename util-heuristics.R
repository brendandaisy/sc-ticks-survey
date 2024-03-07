# --------------------------------------------------------------------------------
# util-heuristics.R---------------------------------------------------------------
# Find designs and compute utility of the two search heuristics used in the-------
# text, selecting visits based on predicted variance and a space-filling one------
# --------------------------------------------------------------------------------
library(INLA)
library(purrr)
library(tidyverse)
library(sf)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

num_visits <- seq(5, 20, 5)
n <- 50

# Variance heuristic--------------------------------------------------------------
# Select sites with largest average variance among species,-----------------------
# only choosing a site one time---------------------------------------------------
pred_mod_sort <- pred_mod |> 
    group_by(date, site) |> 
    # average variance from each species:
    summarise(var_eta=mean(var_eta), .groups="drop") |> 
    arrange(desc(var_eta)) |> 
    group_by(site) |> 
    # only allow each site to be selected once:
    filter(row_number() == 1) |> 
    ungroup()

u_simple_var <- map_dfr(num_visits, ~{
    print(paste0("num. visits =", .x))
    d <- head(pred_mod_sort, .x)
    tibble_row(
        !!!utility(d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod),
        num_loc=.x
    )
})

saveRDS(u_simple_var, "util-results/risk-sd/util-simple-var.rds")

# Space filling heuristic---------------------------------------------------------
# spread out samples over time, with pairs a given dist. away---------------------

#' get_design_spacefill
#'
#' @param sites list of allowed sites
#' @param dates list of allowed dates
#' @param num number of visits
#' @param pairs a list of banned pairs of sites (sites that would be too close)
#'
#' @return tibble of visits in the proposed design
get_design_spacefill <- function(sites, dates, num, pairs) {
    ret <- character(num)
    num_acc <- 1
    ret[num_acc] <- sample(sites, 1)
    sites2 <- sites # for if there are not enough valid locations
    while (num_acc < num) {
        s <- sample(sites, 1)
        banp <- semi_join(tibble(s1=ret, s2=s), pairs, by=c("s1", "s2"))
        if (nrow(banp) > 0) {
            sites <- sites[sites != s]
            if (length(sites) == 0)
                break
            next
        }
        num_acc <- num_acc + 1
        ret[num_acc] <- s
    }
    if (num_acc < num)
        ret[(num_acc+1):num] <- sample(sites2, num-num_acc)
    return(tibble(site=ret, date=choose_dates(dates, length(ret))))
}

# to spread visit times evenly over the 12 months
choose_dates <- function(dates, num) {
    if (length(dates) >= num) {
        return(sample(dates, num, replace=FALSE))
    }
    ret <- sample(dates, replace=FALSE) # shuffle
    c(ret, sample(dates, num-length(dates), replace=TRUE))
}

# get all the park sites in sf format
all_loc_sp <- read_sf("data-proc/parks-design-space.shp") |> 
    distinct(site, geometry)

# for maintaining spatial spread, remove all sites that have already been visited
sp_unvis_loc <- all_loc_sp |> 
    filter(!(site %in% unique(obs_mod$site))) |> 
    distinct(site, geometry)

min_dist <- units::as_units(25000, "m")

# maintain a list of sites that are too close and should not co appear in a design
banned_pairs <- sp_unvis_loc |> 
    cross_join(sp_unvis_loc) |> 
    filter(st_distance(geometry.x, geometry.y, by_element=TRUE) < min_dist) |> 
    transmute(s1=site.x, s2=site.y)

for (r in 2:5) {
    u_sfill_rep <- map_dfr(num_visits, ~{
        print(paste0("num. visits =", .x))
        dsfill <- get_design_spacefill(sp_unvis_loc$site, unique(pred_mod$date), .x, banned_pairs)
        tibble_row(
            !!!utility(dsfill, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod),
            num_loc=.x
        )
    })
    
    saveRDS(u_sfill_rep, paste0("util-results/risk-sd/util-spacefill-", r, ".rds"))
}

## load all spacefill results, pick best and save
u_spacefill <- map_dfr(1:5, ~readRDS(paste0("util-results/risk-sd/util-spacefill-", .x, ".rds")))

u_spacefill <- u_spacefill |> 
    filter(num_loc > 1) |> 
    group_by(num_loc) |> 
    slice_max(risk_sd, n=1) |> 
    ungroup()

saveRDS(u_spacefill, "util-results/risk-sd/util-spacefill.rds")