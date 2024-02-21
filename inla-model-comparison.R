library(INLA)
library(tidyverse)
library(sf)
library(cowplot)
library(ggridges)
# library(pROC)
# library(gridExtra)

source("other-helpers.R")

inla.setOption(inla.mode="classic")

## Format the fixed effects for the INLA formula, with 1=single slope to 4=spline for each species
get_fx <- function(var, level) {
    case_when(
        # level == 0 ~ "",
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

## Combine the fx terms into a formula
get_formula <- function(terms, response="pres") {
    fx <- get_fx(names(terms), terms) |> keep(~str_length(.x) > 1)
    paste0(response, "~ 0 + tick_class +", str_c(fx, collapse="+")) |>
        as.formula()
}

# data to be used in all following analyses.
# covariates are rescaled since no grid locations added in this script
# only land_cover levels appearing in parks are used
parks_obs <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=TRUE)

# Model performance for binary response-------------------------------------------

fx_terms <- rep(1, 7)

names(fx_terms) <- parks_obs |> 
    st_drop_geometry() |> 
    select(land_cover:mean_rh) |> 
    colnames()

# prior choices:
# prec for fixed effect = 0.2 -> an increase in 1 std. dev. has 95% chance of causing a 99% change in risk (each covar has sufficient chance to entirely describe risk)
# prec for random effects AND splines = Gamma(1, 0.1) -> going for more precision than for linear fixed effects, since RWs have so much room to move around per step. A change in "unit" has 2% chance of prec < 0.2
prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.1)))
prec_pri_str <- "list(prec=list(prior='loggamma', param=c(1, 0.1)))" # this one is for get_fx

fx_models <- list(
    slope=get_formula(fx_terms),
    # `slope:species`=get_formula(fx_terms + 1),
    # spline=get_formula(c(land_cover=2, fx_terms + 2)),
    # `spline:species`=get_formula(c(land_cover=2, fx_terms + 3))
    `slope:species`=get_formula(fx_terms + 1),
    spline=get_formula(c(land_cover=2, fx_terms[-1] + 2)),
    `spline:species`=get_formula(c(land_cover=2, fx_terms[-1] + 3))
)

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

mod_comp_res <- expand_grid(fx=fct_inorder(names(fx_models)), rx=fct_inorder(names(rx_models))) |> 
    mutate(
        form=map2(fx, rx, ~update(fx_models[[.x]], rx_models[[.y]])),
        ft=map(form, ~{
            print(.x)
            fit_model(.x, parks_obs, fx_prec=0.2, control_compute=list(mlik=TRUE, dic=TRUE))
        }),
        mlik=map_dbl(ft, ~.x$mlik[1, 1]),
        dic=map_dbl(ft, ~.x$dic$dic),
        cpu=map_dbl(ft, ~.x$cpu.used[4])
    )

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
        # for presentations:
        # strip.text=element_text(size=rel(1.15)),
        # axis.text=element_text(size=rel(1.03)),
        # axis.title=element_text(size=rel(1.15))
    )

ggsave("figs/model-comp-dic.pdf", width=5.2, height=5)

# Best model performance----------------------------------------------------------

best_model <- slice_min(mod_comp_res, dic, n=1)$form[[1]]
best_fit <- fit_model(formula_bed(), parks_obs, 0.2)

## Marginal posteriors for fixed effect coefficients
p1 <- best_fit$marginals.fixed |> 
    summ_fx() |> 
    # mutate(species=str_extract(var, "Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\."), var=fx_labels(var)) |> 
    mutate(var=fx_labels(var)) |> 
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray70', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col="#3a2b52") +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col="#3a2b52") +
    geom_point(size=1.5, col="#9de0c5") +
    labs(x = "Log odds", y="Environmental variable") +
    theme_bw() +
    theme(
        # for presentations:
        # strip.text=element_text(size=rel(1.15)),
        # axis.text=element_text(size=rel(1.03)),
        # axis.title=element_text(size=rel(1.15))
    )

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

ggsave("figs/model-effects.pdf", width=7.5, height=4.7)

as_tibble(best_fit$summary.random$month) |>
    mutate(group=rep(c("Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp."), each=10)) |>
    ggplot(aes(month(ID, label=TRUE, abbr=FALSE), col=group)) +
    geom_line(aes(y=`0.5quant`, group=group)) +
    geom_linerange(aes(ymin=`0.025quant`, ymax=`0.975quant`)) +
    geom_point(aes(y=`0.5quant`)) +
    facet_wrap(~group, nrow=3) +
    labs(x="Month", y="Temporal effect")

ggsave("figs/pres-fixed-effects.pdf", width=5.5, height=7.3)

###
parks_obs <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE)

all_data <- append_pred_data(parks_obs, pred_grid=prep_pred_grid(parks_obs), drop_new_lcc=FALSE)

fit_all <- fit_model(formula_bed(), all_data, fx_prec=0.2)
fit_all$summary.fixed # assert: fixed effects were included in the model correctly?

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

ppg_tmp <- post_pred_grid |> 
    mutate(season=fct_inorder(case_match(as.numeric(month), c(12, 1:2) ~ "Winter", 3:5 ~ "Spring", 6:8 ~ "Summer", 9:11 ~ "Fall"))) |> 
    group_by(site, tick_class, season) |> 
    summarise(mean_pres=mean(mean_pres)) |> 
    ungroup()

ggplot(ppg_tmp, aes(col=mean_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.02) +
    scale_color_viridis_c(option="mako") +
    # new_scale_color() +
    # geom_sf(aes(col=as.factor(pres)), data=pres, size=1.5) +
    # scale_color_manual(values=c(`0`="pink", `1`="red")) +
    facet_grid(tick_class~season, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col=NULL) +
    pred_map_theme +
    theme(
        legend.position="right",
        # strip.text.x=element_text(size=rel(1.05))
        strip.text.y=element_text(size=rel(1.09)),
        legend.title=element_text(size=rel(1.3))
    )

ggsave("astmh-poster/pred-map.pdf", width=5.4, height=3.8)

p2 <- ggplot(post_pred_grid, aes(col=sd_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.02) +
    scale_color_viridis_c(option="magma", breaks=seq(0, 0.4, 0.1), limits=c(0, 0.4)) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Standard deviation") +
    pred_map_theme +
    theme(
        strip.text.x=element_text(size=rel(1.05))
        # strip.text.y=element_text(size=rel(1.05)), 
        # legend.title=element_text(size=rel(1.25))
    )

ggsave("figs/presence-prediction-map.pdf", plot_grid(p1, p2, nrow=2), width=11.5, height=7)

# Reporting the residual variability in collection timing and site----------------

m_sd_mon <- inla.tmarginal(function(x) sqrt(1/x), fit_all$marginals.hyperpar$`Precision for month`)
inla.zmarginal(m_sd_mon)
inla.hpdmarginal(0.9, m_sd_mon)

m_sd_site <- inla.tmarginal(function(x) sqrt(1/x), fit_all$marginals.hyperpar$`Precision for site`)
inla.zmarginal(m_sd_site)
inla.hpdmarginal(0.9, m_sd_site)

# Prior prediction----------------------------------------------------------------
parks_pri_pred <- parks_obs
parks_pri_pred$pres <- NA

ft_pri_pred <- fit_model(formula_bed(), parks_pri_pred, fx_prec=0.2)

summary(ft_pri_pred$summary.linear.predictor$`0.025quant`)

## ROC curves for all observations and each species, based on mean fitted values
# preds <- transmute(parks_model_data, tick_class, pres, pred=best_fit$summary.fitted.values$mean)
# 
# roc <- list(
#     Dermacentor=roc(filter(preds, tick_class == "Dermacentor variabilis"), pres, pred),
#     Ixodes=roc(filter(preds, tick_class == "Ixodes spp."), pres, pred),
#     Amblyomma=roc(filter(preds, tick_class == "Amblyomma americanum"), pres, pred),
#     All=roc(preds, pres, pred)
# )
# 
# names(roc) <- imap_chr(roc, ~paste0(.y, " (AUC=", round(.x$auc, 4), ")"))
# 
# ggroc(roc, alpha=0.75, size=1.3) +
#     theme_bw() +
#     labs(x="Specificity", y="Sensitivity", col=NULL) +
#     theme(
#         legend.position=c(0.7, 0.2), 
#         legend.background=element_rect(color="gray40"), 
#         legend.text=element_text(size=rel(0.6))
#     )
# 
# ggsave("figs/pres-roc.pdf", width=4, height=4)
# 
# ## TODO: Temporal and spatial posterior trends
# as_tibble(fit_pres$summary.random$`jan_min_temp`) |> 
#     mutate(group=rep(1:3, each=102%/%3)) |> 
#     ggplot(aes(ID, `0.5quant`)) +
#     geom_line(aes(col=group, group=group)) +
#     geom_point(data=transmute(pres, ID=jan_min_temp, `0.5quant`=pres))