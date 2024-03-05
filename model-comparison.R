# --------------------------------------------------------------------------------
# model-comparison.R--------------------------------------------------------------
# Implementation of step 1 of BED, where a number of forms for environmental------
# and spatiotemporal effects are compared-----------------------------------------
# --------------------------------------------------------------------------------

library(INLA)
library(tidyverse)
library(sf)
library(cowplot)
# library(ggridges)
# library(gridExtra)

source("other-helpers.R")

# for the `sel` joint posterior feature to work
inla.setOption(inla.mode="classic")

#' get_fx
#' Create different terms for the environmental effects for the INLA formula
#'
#' @param var name(s) of the covariate(s). Should match the column name of the input data
#' @param level format of effect, with 1=single slope to 4=spline for each species
#'
#' @return a string for variable in `var`
get_fx <- function(var, level) {
    case_when(
        level == 0 ~ "",
        level == 1 ~ var,
        level == 2 ~ str_c(var, ":tick_class"),
        level == 3 ~ str_c("f(inla.group(", var, "), model='rw1', hyper=", prec_pri_str, ")"),
        level == 4 ~ str_c(
            "f(inla.group(",
            var,
            "), model='rw1', hyper=",
            prec_pri_str,
            ", group=tcnum, control.group=list(model='iid', hyper=",
            prec_pri_str, "))"
        )
    )
}

#' get_formula
#' Combine the environmental effect terms into an R formula
#'
#' @param terms names of environmental covariates to use
#' @param response name of the response variable, e.g. `pres`
#'
#' @return an INLA formula
get_formula <- function(terms, response="pres") {
    fx <- get_fx(names(terms), terms) |> keep(~str_length(.x) > 1)
    paste0(response, "~ 0 + tick_class +", str_c(fx, collapse="+")) |>
        as.formula()
}

# data to be used in all following analyses.
# covariates are rescaled since no grid locations added in this script
# only land_cover levels appearing in parks are used
parks_obs <- read_parks_sf() |> 
    rescale_covars()

# prior choices:
# prec for fixed effect = 0.2 -> an increase in 1 std. dev. has 95% chance of causing a 99% change in risk (each covar has sufficient chance to entirely describe risk)
# prec for random effects AND splines = Gamma(1, 0.1) -> going for more precision than for linear fixed effects, since RWs have so much room to move around per step
# a change in "unit" has 2% chance of decreasing prec
prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.1)))
prec_pri_str <- "list(prec=list(prior='loggamma', param=c(1, 0.1)))" # this one is for get_fx

# set up different model structures for the environmental effects:
fx_terms <- rep(1, 7)

names(fx_terms) <- parks_obs |> 
    st_drop_geometry() |> 
    select(land_cover:mean_rh) |> 
    colnames()

fx_models <- list(
    slope=get_formula(fx_terms),
    `slope:species`=get_formula(fx_terms + 1),
    spline=get_formula(c(land_cover=2, fx_terms[-1] + 2)), # don't do splines for the land_cover since 
    `spline:species`=get_formula(c(land_cover=2, fx_terms[-1] + 3))
)

# set up different structures for the spatiotemporal effects:
rx_models <- list(
    none=~ .,
    month=~ . + f(month, model="ar1", hyper=prec_pri),
    site=~ . + f(site, model="iid", hyper=prec_pri),
    `month+site`=~. + f(month, model="ar1", hyper=prec_pri) + f(site, model="iid", hyper=prec_pri),
    `month:species`=~ . + f(
        month, model="ar1", hyper=prec_pri,
        group=tcnum, control.group=list(model='iid', hyper=prec_pri)
    ),
    `site:species`=~. + f(
        site, model="iid", hyper=prec_pri,
        group=tcnum, control.group=list(model='iid', hyper=prec_pri)
    ),
    `(month+site):species`=~ . + f(
        month, model="ar1", hyper=prec_pri,
        group=tcnum, control.group=list(model='iid', hyper=prec_pri)
    ) + f(
        site, model="iid", hyper=prec_pri,
        group=tcnum, control.group=list(model='iid', hyper=prec_pri)
    )
)

# combine all combos of the environmental and spatiotempoal effects and compute
# the DIC (and some other metrics)
mod_comp_res <- expand_grid(
    fx=fct_inorder(names(fx_models)), 
    rx=fct_inorder(names(rx_models))
) |> 
    mutate(
        form=map2(fx, rx, ~update(fx_models[[.x]], rx_models[[.y]])),
        ft=map(form, ~{
            print(.x)
            fit_model(.x, parks_obs, fx_prec=0.2, control_compute=list(mlik=TRUE, dic=TRUE))
        }),
        mlik=map_dbl(ft, ~.x$mlik[1, 1]), # marginal likelihood. not great for comparing these types of models
        dic=map_dbl(ft, ~.x$dic$dic), # deviance information criterion
        cpu=map_dbl(ft, ~.x$cpu.used[4]) # compute time needed to fit each model
    )

# Figure S1-----------------------------------------------------------------------
# Plot comparing the permance of the 28 candidate models--------------------------
ggplot(mod_comp_res, aes(fx, rx)) +
    geom_tile(aes(fill=dic), col="white") +
    geom_text(aes(label=round(dic, 1)), col="gray80", size=4.9) +
    geom_tile(fill=NA, col="orange", data=slice_min(mod_comp_res, dic, n=1), linewidth=1) +
    # facet_grid(.~jsdm, scales="free_x", space="free", labeller=labeller(jsdm=label_wrap_gen(15))) +
    # scale_fill_viridis_c(option="viridis", values=scales::rescale(sort(mod_comp_res$dic))^2) +
    scale_fill_viridis_c(option="mako") +
    labs(x="Environmental effects", y="Spatiotemporal effects", fill="DIC") +
    theme_bw() +
    scale_x_discrete(expand=expansion(add=c(0.55, 0))) +
    scale_y_discrete(expand=expansion(add=c(0.55, 0.51))) +
    theme(
        panel.border=element_blank(),
        axis.ticks=element_line(color="gray90")
    )

# Plots analyzing results from the best model-------------------------------------
best_model <- slice_min(mod_comp_res, dic, n=1)$form[[1]]
best_fit <- fit_model(formula_bed(), parks_obs, 0.2)

# Figure 2------------------------------------------------------------------------
# Marginal posteriors for fixed effect coefficients
p1 <- best_fit$marginals.fixed |> 
    summ_fx() |> 
    mutate(var=fx_labels(var)) |> 
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray70', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col="#3a2b52") +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col="#3a2b52") +
    geom_point(size=1.5, col="#9de0c5") +
    labs(x = "Log odds", y="Environmental variable") +
    theme_bw()

# Temporal posterior trends
temp_post_samp <- tibble(
    marg=best_fit$marginals.random$month,
    month=rep(lubridate::month(3:12, label=TRUE, abbr=TRUE), 3),
    group=rep(c("Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp."), each=10)
) |> 
    mutate(rmarg=map(marg, ~inla.rmarginal(1000, .x)))

temp_post_samp <- temp_post_samp |> 
    mutate(mean=map_dbl(rmarg, mean))

p2 <- ggplot(unnest(temp_post_samp, rmarg), aes(rmarg, month, fill=group)) +
    geom_density_ridges(scale=1, alpha=0.35, col=NA) +
    geom_line(aes(x=mean, group=group, col=group), data=temp_post_samp, orientation="y", linetype="dotdash", linewidth=1) +
    geom_point(aes(x=mean, col=group), data=temp_post_samp, size=1.1) +
    theme_bw() +
    xlim(-5, 5) +
    coord_flip() +
    labs(y="Month", x="Temporal effect") +
    scale_fill_manual(values=c("#3888a6", "#3a2b52", "#9de0c5")) +
    scale_color_manual(values=c("#3888a6", "#3a2b52", "#9de0c5")) +
    theme(legend.position="bottom", legend.direction="vertical", legend.title=element_blank())

plot_grid(p1, p2, nrow=1, rel_heights=c(1, 0.86), rel_widths=c(1.1, 1), labels=c("A", "B"))

###
pred_grid <- prep_pred_grid(parks_obs)

all_data <- read_parks_sf() |> # original parks data without scaled covariates
    append_pred_data(pred_grid=pred_grid, drop_new_lcc=FALSE)

fit_all <- fit_model(formula_bed(), all_data, fx_prec=0.2)

post_pred_grid <- all_data |>
    mutate(
        month=lubridate::month(month, label=TRUE, abbr=FALSE), # month to english
        tick_class=str_replace(tick_class, " ", "\n"),
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd
    ) |>
    filter(is.na(pres) & str_detect(site, "g\\d+"))

pred_map_theme <- theme_bw() +
    theme(
        axis.text=element_blank(), axis.ticks=element_blank(),
        panel.spacing=unit(0, "mm"),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
        legend.position="bottom",
        legend.box.margin=unit(c(0, 0, 0, 0), "mm")
    )

ggplot(post_pred_grid, aes(col=mean_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.02) +
    scale_color_viridis_c(option="mako") +
    # scale_color_viridis_c(option="magma", breaks=seq(0, 0.4, 0.1), limits=c(0, 0.4)) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Probability tick present") +
    pred_map_theme
