library(INLA)
library(tidyverse)
library(sf)
library(lubridate)
# library(MASS, exclude=c("select"))
# library(furrr)
# library(gridExtra)
# library(DiceOptim)
# library(fields)
library(caret)
library(ggfortify)
library(cowplot)

source("other-helpers.R")
source("utility-helpers.R")

inla.setOption(inla.mode="classic")

# Data preparation----------------------------------------------------------------
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

num_loc <- c(1, seq(5, 20, 5))
n <- 50

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

u_simple_var <- map_dfr(num_loc, ~{
    print(paste0("num_loc=", .x))
    d <- head(pred_mod_sort, .x)
    tibble_row(
        !!!utility(d, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod),
        num_loc=.x
    )
})

saveRDS(u_simple_var, "util-results/risk-sd/util-simple-var.rds")

# Space-filling strategy----------------------------------------------------------
# spread out samples over time, with pairs a given dist. away---------------------
choose_dates <- function(dates, num) {
    if (length(dates) >= num) {
        return(sample(dates, num, replace=FALSE))
    }
    ret <- sample(dates, replace=FALSE) # shuffle
    c(ret, sample(dates, num-length(dates), replace=TRUE))
}

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

sp_unvis_loc <- all_loc_sp |> 
    filter(!(site %in% unique(obs_mod$site))) |> 
    distinct(site, geometry)

min_dist <- units::as_units(25000, "m")

banned_pairs <- sp_unvis_loc |> 
    cross_join(sp_unvis_loc) |> 
    filter(st_distance(geometry.x, geometry.y, by_element=TRUE) < min_dist) |> 
    transmute(s1=site.x, s2=site.y)

for (r in 2:5) {
    u_sfill_rep <- map_dfr(num_loc[-1], ~{
        print(paste0("num_loc=", .x))
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

# Local utility stuff-------------------------------------------------------------

local_utility <- function(d, known_df, pred_df=NULL, copy=FALSE) {
# TODO: revisit at some point whether known_df can be removed, or can be used for only one tick species
# probably can, but confusing especially with the dummy matrix
    d_df <- if (!is.null(pred_df)) inner_join(pred_df, d, by=c("date", "site"), copy=copy) else d
    d_df <- bind_rows(d_df, known_df)
    
    X <- dummyVars(
        pres ~ -1 + tick_class + tick_class:tree_canopy + tick_class:elevation +
            tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh,
        data=d_df,
        fullRank=TRUE
    )
    X <- predict(X, d_df)
    
    I <- t(X) %*% diag(d_df$mean_pres * (1 - d_df$mean_pres)) %*% X
    return(log(det(I))) # want to MAXIMIZE information to minimize uncertainty
}

coarsen_x <- function(x, length=100) {
    xc <- seq(min(x), max(x), length.out=length)
    int <- findInterval(x, xc)
    xc[int]
}

get_matching_x <- function(pred_df, x) {
    ret <- inner_join(pred_df, x, by=colnames(x))
    if (nrow(ret) <= 1)
        return(ret)
    
    ret |> 
        select(date, site) |> 
        slice_sample(n=1)
}

opt_coord_exchange <- function(known_df, pred_df, d_so_far=tibble(), n_iter=1, num_loc=1, x_res=100) {
    
    pred_df_coarse <- mutate(pred_df, across(tree_canopy:mean_rh, coarsen_x, length=x_res))
    
    d_cur <- pred_df_coarse |> 
        select(date, site) |> 
        slice_sample(n=num_loc)
    
    x_choices <- pred_df_coarse |> 
        select(land_cover:mean_rh) |> 
        map(unique)
    
    x_names <- names(x_choices)
    d_df_init <- inner_join(pred_df_coarse, d_cur, by=c("date", "site"))
    x_cur <- d_df_init[1:num_loc, x_names]
    
    u_best <- utility(d_df_init, known_df, n=10, u_only=TRUE)
    
    for (iter in 1:n_iter) {
        for (s in 1:num_loc) { # modify each location in the current batch at a time
            loc_best <- d_cur[s,]
            for (xi in seq_along(x_choices)) { # modify each coordinate at a time
                x_best <- x_cur[,xi]
                nmatch <- 0
                for (x_val in x_choices[[xi]]) {
                    x_cur[,xi] <- x_val
                    loc_new <- get_matching_x(pred_df_coarse, x_cur)
                    if (nrow(loc_new) == 0)
                        next
                    nmatch <- nmatch + 1
                    d_cur[s,] <- loc_new
                    u <- utility(d_cur, known_df, full_df=pred_df_coarse, n=10, u_only=TRUE)
                    if (u > u_best) {
                        u_best <- u
                        x_best <- x_val
                        loc_best <- loc_new
                        print(paste0("Success! U=", u))
                        print(paste0(names(x_choices)[xi], ": ", x_best))
                        print(loc_best)
                    }
                }
                x_cur[,xi] <- x_best
                print(paste0("Num. valid: ", nmatch))
            }
            d_cur[s,] <- loc_best
            print(d_cur)
        }
    }
    return(d_cur)
}

d_xco <- opt_coord_exchange(parks_mod, grid_mod, n_iter=5, num_loc=1, x_res=10)

u_local_random <- map_dfr(rep(num_loc, each=100), ~{
    d <- grid_mod |> 
        select(date, site) |> 
        slice_sample(n=.x)
    
    tibble_row(size=.x, d=list(d), u=local_utility(d, parks_mod, grid_mod, copy=TRUE))
})

ggplot(u_local_random, aes(size, u)) +
    geom_point() +
    annotate("point", x=5, y=local_utility(d_xco, parks_mod, grid_mod), col="red")

d <- grid_mod |> 
    select(date, site) |> 
    slice_sample(n=1)

local_utility(d, parks_mod, grid_mod)
local_utility(tibble(), parks_mod)

tmp <- map_dbl(1:100, ~local_utility(d, parks_mod, grid_mod))

u_xco <- utility(d_xco, parks_mod, full_df=grid_mod, n=50)

tmp <- mutate(grid_mod, across(tree_canopy:mean_rh, coarsen_x, length=10))

distinct(tmp, across(land_cover:mean_rh))

# "One-off" design strategies-----------------------------------------------------

## Visit all 30 parks again in December/June

d_parks_dec <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-12-01"))
    
u_parks_dec <- utility(d_parks_dec, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)

d_parks_jun <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-06-01"))

u_parks_jun <- utility(d_parks_jun, obs_mod, n=n, pred_df=pred_mod, util_fun=util_risk_sd, risk_df=risk_mod)

## Reuse the same schedule from 2021
d_resamp <- filter(obs_mod, lubridate::year(date) == 2021)
d_resamp$pres <- NA
u_resamp <- utility(d_resamp, obs_mod, n=n, util_fun=util_risk_sd, risk_df=risk_mod)

u_one_offs <- tibble(
    strat=c("    Revisit in\n    June (30)", "    Revisit in\n    December (30)", "    Repeat 2021\n    schedule (111)"), 
    design=list(distinct(d_parks_jun, date, site), distinct(d_parks_dec, date, site), distinct(d_resamp, date, site)),
    utility=c(u_parks_jun$risk_sd, u_parks_dec$risk_sd, u_resamp$risk_sd)
)

ofs_dopt <- readRDS("util-results/util-one-offs.rds")
ofs_rsd <- readRDS("util-results/risk-sd/util-one-offs.rds")

# Plot the results comparing all strategies---------------------------------------

strats <- c("util-random", "util-simple-var", "util-spacefill", "util-sim-ann", "util-coord-ex")
names(strats) <- c("Random", "Variance", "Space-filling", "Simulated Annealing", "Exchange")

res_dopt <- strats |> 
    map_dfr(~mutate(readRDS(paste0("util-results/", .x, ".rds")), strat=.x)) |> 
    rename(utility=dopt) |> 
    filter(num_loc > 1) |> 
    mutate(strat=fct_relevel(fct_recode(strat, !!!strats), "Variance", after=Inf))

res_dopt_random <- filter(res_dopt, strat == "Random")

res_rsd <- strats |> 
    map_dfr(~mutate(readRDS(paste0("util-results/risk-sd/", .x, ".rds")), strat=.x)) |> 
    rename(utility=risk_sd) |> 
    filter(num_loc > 1) |> 
    mutate(strat=fct_relevel(fct_recode(strat, !!!strats), "Variance", after=Inf))

res_rsd_random <- filter(res_rsd, strat == "Random")

hline_col <- c("#3798a9", "#7fd7b8", "#98d180")
    
gg_dopt <- ggplot(res_dopt_random, aes(factor(num_loc), utility)) +
    geom_hline(yintercept=ofs_dopt$utility, col=hline_col, alpha=0.65, linewidth=1, linetype="dashed") +
    geom_violin(fill="#9ac9e7", col="gray70", alpha=0.75)  +
    geom_point(size=0.84, col="#9ac9e7") +
    coord_cartesian(clip="off") +
    labs(x=NULL, y=NULL, col="Search strategy") +
    theme_bw() +
    theme(plot.margin=unit(c(t=0.1, r=0.1, b=0.1, l=0), "in"))

gg_rsd <- ggplot(res_rsd_random, aes(factor(num_loc), utility)) +
    geom_hline(yintercept=ofs_rsd$utility, col=hline_col, alpha=0.65, linewidth=1, linetype="dashed") +
    geom_violin(fill="#9ac9e7", col="gray70", alpha=0.75)  +
    geom_point(size=0.84, col="#9ac9e7") +
    coord_cartesian(clip="off") +
    annotate(
        "text",
        label=ofs_rsd$strat,
        # y=ofs_rsd$utility + c(-0.06, 0.06, 0),
        y=ofs_rsd$utility + c(0.002, -0.002, 0),
        col=hline_col, x=4.45, hjust="left", size=3
    ) +
    labs(x=NULL, y=NULL, col=NULL) +
    theme_bw() +
    theme(plot.margin=unit(c(t=0.1, r=1, b=0.1, l=0), "in"))


plot_searches <- function(gg, u_res, vars) {
    u_sub <- filter(u_res, strat != "Random", strat %in% vars) |> 
        mutate(strat=fct_rev(strat))
    
    gg +
        geom_line(aes(col=strat, group=strat), u_sub) +
        geom_point(aes(col=strat, shape=strat, fill=strat), u_sub, size=1.1, show.legend=FALSE) +
        scale_color_manual(values=c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc")[seq_along(vars)]) +
        scale_fill_manual(values=c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc")[seq_along(vars)]) +
        scale_shape_manual(values=21:24, guide="none") +
        guides(color=guide_legend(title.position="top")) +
        theme(legend.position="none")
}

gg_dopt <- plot_searches(gg_dopt, res_dopt, names(strats))
gg_risk_sd <- plot_searches(gg_rsd, res_rsd, names(strats))
legend <- get_legend(gg_dopt + theme(legend.position="top", legend.justification="left"))

plot_grid(
    legend,
    # plot_grid(gg_dopt, NULL, rel_widths=c(1, 0.3)),
    plot_grid(gg_dopt, NULL, gg_risk_sd, rel_widths=c(1, 0.05, 1.39), nrow=1),
    # plot_grid(NULL, legend, rel_widths=c(0.2, 1)),
    nrow=2, rel_heights=c(0.15, 1)
)

ggsave("figs/figure 4/util-results-a.pdf", width=6.2, height=3.8)

# Visualizing designs using dimension reduction-----------------------------------
#  Figure S3 of TTBDis paper------------------------------------------------------
library(FactoMineR)
library(factoextra)

famd_searches <- function(res, covars, num_loc=10) {
    covar_search <- res |> 
        filter(strat != "Random", !!num_loc == num_loc) |> 
        unnest(cols=c(design)) |> 
        mutate(month=ifelse(is.na(month), lubridate::month(date), month)) |> 
        left_join(covars, by=c("site", "month"))
    
    pcres <- FAMD(select(covars, land_cover:mean_rh), graph=FALSE)
    list(
        res=pcres,
        obs=predict(pcres, obs_mod)$coord |> as_tibble() |> mutate(strat="Initial Data"),
        search=predict(pcres, covar_search)$coord |> as_tibble() |> mutate(strat=covar_search$strat)
    )
}

plot_famd_searches <- function(famd, d1, d2, axes=c(1, 2)) {
    fviz_famd_ind(famd$res, axes=axes, label="none", col.quali.var="white", col.ind="gray60", alpha.ind=0.7, shape.ind=1, pointsize=0.75) +
        geom_point(aes({{d1}}, {{d2}}), data=famd$obs, col="gray60", alpha=0.7, size=0.75) +
        geom_point(aes({{d1}}, {{d2}}, fill=fct_rev(strat), shape=fct_rev(strat)), data=famd$search, size=1.7, col="white") +
        scale_fill_manual(values=c("#f5b43d", "#ff856d", "#c63df5", "#f53dbc"), guide="none") +
        scale_shape_manual(values=21:24, guide="none") +
        coord_cartesian(clip="off") +
        labs(title=NULL) +
        theme_bw() +
        theme(plot.margin=unit(c(t=0, r=0.1, b=0.1, l=0.1), "in"))
}

covars <- pred_mod |> 
    select(date, month, site, land_cover:mean_rh) |> 
    distinct()

famd1 <- famd_searches(res_dopt, covars)
famd2 <- famd_searches(res_rsd, covars)

famd_grid <- plot_grid(
    plot_famd_searches(famd1, `Dim 1`, `Dim 2`) + labs(title="1st design criterion"),
    plot_famd_searches(famd2, `Dim 1`, `Dim 2`) + labs(title="2nd design criterion"),
    plot_famd_searches(famd1, `Dim 3`, `Dim 4`, c(3, 4)),
    plot_famd_searches(famd2, `Dim 3`, `Dim 4`, c(3, 4)),
    nrow=2, labels="AUTO"
)

plot_grid(legend, famd_grid, nrow=2, rel_heights=c(0.15, 1))

ggsave("figs/supp/util-res-famd.pdf", width=5.75, height=6.2)

# Exploring other aspects of repeated sites, visits, etc--------------------------
dopt_unnest <- res_dopt |>
    filter(!(strat %in% c("Random", "Variance")), num_loc == 20) |> 
    unnest(cols=c(design)) |> 
    # mutate(month=lubridate::month(ifelse(is.na(month), lubridate::month(date), month), label=TRUE)) |> 
    mutate(month=ifelse(is.na(month), lubridate::month(date), month)) |> 
    left_join(pred_mod, by=c("site", "month"), multiple="first") |> 
    mutate(month=month(month, label=TRUE, abbr=FALSE))

rsd_unnest <- res_rsd |>
    filter(!(strat %in% c("Random", "Variance")), num_loc == 20) |> 
    unnest(cols=c(design)) |> 
    # mutate(month=lubridate::month(ifelse(is.na(month), lubridate::month(date), month), label=TRUE)) |> 
    mutate(month=ifelse(is.na(month), lubridate::month(date), month)) |> 
    left_join(pred_mod, by=c("site", "month"), multiple="first") |> 
    mutate(month=month(month, label=TRUE, abbr=FALSE))

# save a CSV of the best search strategies
all_loc_sp <- read_sf("data-proc/parks-design-space.shp") |> 
    distinct(site, geometry)

best_res <- bind_rows(
    mutate(slice_max(dopt_unnest, utility), criterion="first"),
    mutate(slice_max(rsd_unnest, utility), criterion="second")
) |> 
    select(site, month, num_loc, strat, criterion, utility)

inner_join(all_loc_sp, best_res, multiple="all") |> 
    arrange(criterion, month) %>%
    mutate(long=st_coordinates(.)[,1], lat=st_coordinates(.)[,2]) |> 
    st_drop_geometry() |> 
    write_csv("ttbd-submission/data-files/best-designs.csv")

# most sampled land cover classes
# risk_mod |> 
#     count(land_cover, sort=TRUE) |> 
#     mutate(n=n/sum(n))
# 
# ggplot(dopt_unnest, aes(month, fill=fct_rev(strat))) +
#     # geom_point(position=position_dodge2(width=0.5), size=1.3) +
#     geom_bar(position=position_dodge2(preserve="single", padding=0.2)) +
#     scale_fill_manual(values=c("#ff856d", "#c63df5", "#f53dbc"), guide="none") +
#     scale_y_continuous(breaks=0:5) +
#     scale_x_discrete(guide=guide_axis(angle=45)) +
#     labs(x=NULL, y="Count") +
#     theme_bw() +
#     theme(panel.grid.minor=element_blank())
# 
# ggsave("figs/figure 4/timepoints-dopt.pdf", width=3, height=1.4)

###


# dopt_best <- slice_max(res_dopt, utility, n=1) |> 
#     unnest(design) |> 
#     mutate(month=lubridate::month(month, label=TRUE, abbr=FALSE))
# 
# rsd_best <- slice_max(res_rsd, utility, n=1) |> 
#     unnest(design) |> 
#     mutate(month=lubridate::month(date, label=TRUE, abbr=FALSE))

dopt_sp <- inner_join(all_loc_sp, dopt_unnest, multiple="all", copy=TRUE)
    # mutate(strat=factor(strat, labels=c("1st design criterion", "2nd design criterion"))) |> 
    # st_jitter(factor=0.02)

rsd_sp <- inner_join(all_loc_sp, rsd_unnest, multiple="all", copy=TRUE)
    # mutate(strat=factor(strat, labels=c("1st design criterion", "2nd design criterion"))) |> 
    # st_jitter(factor=0.02)

sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(st_crs(parks_obs)) |> 
    st_union() |> 
    nngeo::st_remove_holes()

design_map_theme <- theme_bw() +
    theme(
        axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
        panel.spacing=unit(0.55, "mm"),
        # strip.text=element_text(size=rel(0.5)),
        # strip.background=element_rect(fill="gray70", linewidth=0),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
        panel.grid=element_blank(),
        legend.position="none"
    )

p2 <- ggplot() +
    geom_sf(data=sc_state, fill=NA, col="gray70", linewidth=0.7, alpha=0.9) +
    geom_sf(data=init_visits, col="gray70", alpha=0.8, size=1.2) +
    geom_sf(aes(fill=strat, shape=strat), rsd_sp, size=2.1, col="white", alpha=0.7) +
    # geom_sf_text(aes(label=site), col=col, alpha=0.8, size=1.6, nudge_y=-0.2, nudge_x=0.2) +
    facet_wrap(~month, drop=FALSE, nrow=2) +
    scale_fill_manual(values=c("#ff856d", "#c63df5", "#f53dbc")) +
    scale_shape_manual(values=22:24, guide="none") +
    design_map_theme

plot_grid(p1, NULL, p2, nrow=3, rel_heights=c(1, 0.09, 1))

ggsave("figs/figure 4/designs-map.pdf", width=6, height=4.5)

## "how many proposed design points intersected with the initial collections?"
init_visits |> 
    distinct(site, month) |>  # get rid of any repeated visits in the same month
    inner_join(best_res) |> 
    select(site, month, strat, criterion)

## "how many design points intersected between strategies?"
rsd_unnest |> 
    distinct(site, month)
