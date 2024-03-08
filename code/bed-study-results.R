# --------------------------------------------------------------------------------
# bed-study-results.R-------------------------------------------------------------
# analyze and plot the results from each of the different search techniques-------
# --------------------------------------------------------------------------------
library(INLA)
library(tidyverse)
library(sf)
library(lubridate)
library(ggfortify)
library(cowplot)
library(FactoMineR)
library(factoextra)

source("other-helpers.R")
source("utility-helpers.R")

# Data prep. Only the dfs for the second (risk-sd) criterion should be used
df_list <- readRDS("data-proc/bo-risk-sd-dfs.rds")
obs_mod <- df_list$obs
pred_mod <- df_list$pred
risk_mod <- df_list$risk_grid

# Plot the results comparing all strategies---------------------------------------

# combine and save results from the different search techniques, or load them:

# strats <- c("util-random", "util-simple-var", "util-spacefill", "util-sim-ann", "util-coord-ex")
# names(strats) <- c("Random", "Variance", "Space-filling", "Simulated Annealing", "Exchange")
# 
# res_dopt <- strats |> 
#     map_dfr(~mutate(readRDS(paste0("util-results/", .x, ".rds")), strat=.x)) |> 
#     rename(utility=dopt) |> 
#     filter(num_loc > 1) |> 
#     mutate(strat=fct_relevel(fct_recode(strat, !!!strats), "Variance", after=Inf))
# 
# res_rsd <- strats |> 
#     map_dfr(~mutate(readRDS(paste0("util-results/risk-sd/", .x, ".rds")), strat=.x)) |> 
#     rename(utility=risk_sd) |> 
#     filter(num_loc > 1) |> 
#     mutate(strat=fct_relevel(fct_recode(strat, !!!strats), "Variance", after=Inf))
# 
# saveRDS(res_dopt, "util-results/bed-results-dopt.rds")
# saveRDS(res_rsd, "util-results/bed-results-risk-sd.rds")

res_dopt <- readRDS("util-results/bed-results-dopt.rds")
res_rsd <- readRDS("util-results/bed-results-risk-sd.rds")

res_dopt_random <- filter(res_dopt, strat == "Random")
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
    plot_grid(gg_dopt, NULL, gg_risk_sd, rel_widths=c(1, 0.05, 1.39), nrow=1),
    nrow=2, rel_heights=c(0.15, 1)
)

# Visualizing designs using dimension reduction (Fig S3)--------------------------
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
