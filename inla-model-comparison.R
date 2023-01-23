library(INLA)
library(tidyverse)
library(sf)
library(pROC)
library(gridExtra)

source("other-helpers.R")

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
parks_data <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=TRUE)

# Model performance for binary response-------------------------------------------

fx_terms <- rep(1, 7)

names(fx_terms) <- parks_data |> 
    st_drop_geometry() |> 
    select(land_cover:mean_rh) |> 
    colnames()

# TODO: detailed writeup on these choices
# a bit more precision than for linear fixed effects, since with scaled covars less chance to "move around" there
prec_pri <- list(prec=list(prior="loggamma", param=c(2, 0.1)))
prec_pri_str <- "list(prec=list(prior='loggamma', param=c(2, 0.1)))" # this one is for get_fx

fx_models <- list(
    slope=get_formula(fx_terms),
    `slope:species`=get_formula(fx_terms + 1),
    spline=get_formula(c(land_cover=2, fx_terms + 2)),
    `spline:species`=get_formula(c(land_cover=2, fx_terms + 3))
    # `slope:species`=get_formula(fx_terms + 1),
    # spline=get_formula(c(land_cover=2, fx_terms[-1] + 2)),
    # `spline:species`=get_formula(c(land_cover=2, fx_terms[-1] + 3))
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

# TODO: write this (and everything) up: I also repeated this analysis without land_cover,
# cause it seemed like the model was almost singular due to large cov_mat in those directions,
# but the models were a bit worse (even by ML)

mod_comp1 <- expand_grid(fx=fct_inorder(names(fx_models)), rx=fct_inorder(names(rx_models))) |> 
    mutate(
        form=map2(fx, rx, ~update(fx_models[[.x]], rx_models[[.y]])),
        ft=map(form, ~{
            print(.x)
            fit_model(.x, parks_data, fx_prec=0.2, control_compute=list(mlik=TRUE, dic=TRUE))
        }),
        mlik=map_dbl(ft, ~.x$mlik[1, 1]), dic=map_dbl(ft, ~.x$dic$dic)
    )

mod_comp2 <- expand_grid(fx=fct_inorder(names(fx_models)[1:2]), rx=fct_inorder(names(rx_models))) |> 
    mutate(
        form=map2(fx, rx, ~update(
            update(fx_models[[.x]], ~. + f(
                id, model="iid3d", n=nrow(parks_data), 
                hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3))))
            )),
            rx_models[[.y]]
        )),
        ft=map(form, ~{
            print(.x)
            fit_model(.x, parks_data, fx_prec=0.2, control_compute=list(mlik=TRUE, dic=TRUE))
        }),
        mlik=map_dbl(ft, ~.x$mlik[1, 1]), dic=map_dbl(ft, ~.x$dic$dic)
    )

mod_comp_res <- bind_rows(
    mutate(mod_comp1, jsdm="independent residuals"), mutate(mod_comp2, jsdm="joint species residuals")
) |> 
    mutate(cpu=map_dbl(ft, ~.x$cpu.used[4]))

ggplot(mod_comp_res, aes(fx, rx)) +
    geom_tile(aes(fill=mlik), col="white") +
    # geom_tile(fill=NA, col="orange", data=slice_max(mod_comp_res, mlik, n=1), size=1.3) +
    geom_text(aes(label=round(mlik, 1)), col="gray80", size=4.9) +
    facet_grid(.~jsdm, scales="free_x", space="free", labeller=labeller(jsdm=label_wrap_gen(15))) +
    scale_fill_viridis_c(option="viridis", values=scales::rescale(sort(mod_comp_res$mlik))^2) +
    # scale_fill_viridis_c(option="mako") +
    labs(x="Fixed effects", y="Random effects", fill="Marginal Likelihood") +
    theme_bw() +
    scale_x_discrete(expand=expansion()) +
    scale_y_discrete(expand=expansion()) +
    theme(
        panel.spacing=unit(0,"lines"), 
        # for presentations:
        strip.text=element_text(size=rel(1.15)),
        axis.text=element_text(size=rel(1.03)),
        axis.title=element_text(size=rel(1.15))
    )

ggsave("bdhsc-conf-pres/model-comp-mlik.pdf", width=7.4, height=5.2)

# Best model performance----------------------------------------------------------

best_model <- slice_min(mod_comp_res, dic, n=1)$form[[1]]
best_fit <- fit_model(best_model, parks_data, 0.2)

## Marginal posteriors for fixed effect coefficients
best_fit$marginals.fixed |> 
    summ_fx() |> 
    mutate(species=str_extract(var, "Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\."), var=fx_labels(var)) |> 
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col="#3888a6") +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col="#3888a6") +
    geom_point(size = 1.5, col="#dd4a2c") +
    facet_wrap(~species, nrow=3) +
    labs(x = 'Log odds', y = NULL) +
    theme_bw() +
    theme(
        # for presentations:
        strip.text=element_text(size=rel(1.15)),
        axis.text=element_text(size=rel(1.03)),
        axis.title=element_text(size=rel(1.15))
    )

ggsave("bdhsc-conf-pres/pres-fixed-effects.pdf", width=5.5, height=7.3)

###

all_data <- read_parks_sf(drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) |> 
    append_pred_grid("geo-files/covar-grid-16km.shp")

fit_all <- fit_model(formula_jsdm(all_data, rx_models$`(month+site):species`), all_data, fx_prec=0.2)
fit_all$summary.fixed # assert: fixed effects were included in the model correctly?

post_pred_grid <- all_data |>
    mutate(
        month=lubridate::month(month, label=TRUE), # month to english
        tick_class=str_replace(tick_class, " ", "\n"),
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd,
        var_eta=fit_all$summary.linear.predictor$sd^2
    ) |> 
    filter(is.na(pres))

pred_map_theme <- theme_bw() +
    theme(
        axis.text=element_blank(), axis.ticks=element_blank(),
        panel.spacing=unit(0, "mm"),
        plot.margin=unit(c(0, 0, 0, 0), "mm"),
        legend.position="bottom",
        legend.box.margin=unit(c(0, 0, 0, 0), "mm")
    )

p1 <- ggplot(post_pred_grid, aes(col=mean_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c(option="mako") +
    # new_scale_color() +
    # geom_sf(aes(col=as.factor(pres)), data=pres, size=1.5) +
    # scale_color_manual(values=c(`0`="pink", `1`="red")) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Probability tick present") +
    pred_map_theme +
    theme(
        strip.text.x=element_text(size=rel(1.2)), 
        strip.text.y=element_text(size=rel(1.05)), 
        legend.title=element_text(size=rel(1.25))
    )

p2 <- ggplot(post_pred_grid, aes(col=sd_pres)) +
    geom_sf(shape=15, alpha=0.95, size=1.1) +
    scale_color_viridis_c(option="magma", breaks=seq(0, 0.4, 0.1), limits=c(0, 0.4)) +
    facet_grid(tick_class~month, labeller=as_labeller(c(tick_class=label_wrap_gen))) +
    labs(col="Standard deviation") +
    pred_map_theme +
    theme(
        strip.text.x=element_text(size=rel(1.2)), 
        strip.text.y=element_text(size=rel(1.05)), 
        legend.title=element_text(size=rel(1.25))
    )

ggsave("bdhsc-conf-pres/presence-prediction-map.pdf", grid.arrange(p1, p2, nrow=2), width=12, height=8)

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