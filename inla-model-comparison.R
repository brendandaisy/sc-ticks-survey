library(INLA)
library(tidyverse)
library(sf)
library(pROC)

source("other-helpers.R")

## Format the fixed effects for the INLA formula, with 1=single slope to 4=spline for each species
get_fx <- function(var, level) {
    case_when(
        # level == 0 ~ "",
        level == 1 ~ var,
        level == 2 ~ str_c(var, ":tick_class"),
        level == 3 ~ str_c("f(inla.group(", var, "), model='rw1', hyper=", 2, ")"),
        level == 4 ~ str_c(
            "f(inla.group(",
            var,
            "), model='rw1', hyper=",
            2,
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

# Initial data preparation--------------------------------------------------------
# in all analyses, A. maculatum, all larvae, and the min_temp---------------------
# covariate are removed-----------------------------------------------------------

#TODO: with the new covariates, remove jan_min_temp as well!!!
#TODO: or not? Now seems like not nec. with the == larva bug removed (8/16)

parks <- st_read("geo-files/parks-with-covars.shp") |> 
    rename(
        life_stage=lif_stg, land_cover=lnd_cvr, tree_canopy=tr_cnpy, elevation=elevatn, 
        min_temp=min_tmp, max_temp=max_tmp, precipitation=prcpttn, jan_min_temp=jn_mn_t
    ) |> 
    # Larvae must absolute be filtered since they were always assigned A. americanum!!!
    filter(species != "maculatum", life_stage != "larva") |> 
    select(-min_temp)

#TODO: rerun this to double check and make some changes: 1. give land_cover proper labels

# rescale covariates and format a few things
parks_rescale <- parks |> 
    mutate(
        across(tree_canopy:mean_rh, ~(.x - mean(.x)) / sd(.x)), # center and scale
        land_cover=land_cover_labels(land_cover) |> fct_drop(),
        tick_class=ifelse(genus == "Ixodes", str_c(genus, " spp."), str_c(genus, " ", species))
    )

# data to be used in all following analyses:
parks_model_data <- parks_rescale |> 
    group_by(across(c(date, month, site, tick_class, land_cover:mean_rh))) |>
    summarize(pres=ifelse(sum(count) > 0, 1L, 0L), abun=sum(count), .groups="drop") |> 
    mutate(
        id=1:n(), 
        tcnum=case_when(
            tick_class == "Amblyomma americanum" ~ 1L, 
            tick_class == "Dermacentor variabilis" ~ 2L, 
            TRUE ~ 3L
        )
    ) |> 
    arrange(tcnum)

# Model performance for binary response-------------------------------------------

fx_terms <- rep(1, 7)

names(fx_terms) <- parks_model_data |> 
    st_drop_geometry() |> 
    select(land_cover:mean_rh) |> 
    colnames()

prec_pri <- list(prec=list(prior="loggamma", param=c(1, 0.5)))
prec_pri_str <- "list(prec=list(prior='loggamma', param=c(1, 0.5)))" # this one is for get_fx

fx_models <- list(
    slope=get_formula(fx_terms),
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

mod_comp1 <- expand_grid(fx=fct_inorder(names(fx_models)), rx=fct_inorder(names(rx_models))) |> 
    mutate(
        form=map2(fx, rx, ~update(fx_models[[.x]], rx_models[[.y]])),
        mlik=map_dbl(form, ~{
            print(.x)
            fit_model(.x, data=parks_model_data, response="pres")$mlik[1, 1]
        })
    )

mod_comp2 <- expand_grid(fx=fct_inorder(names(fx_models)[1:2]), rx=fct_inorder(names(rx_models))) |> 
    mutate(
        form=map2(fx, rx, ~update(
            update(fx_models[[.x]], ~. + f(
                id, model="iid3d", n=nrow(pres), 
                hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3))))
            )),
            rx_models[[.y]]
        )),
        mlik=map_dbl(form, ~fit_model(.x, data=parks_model_data, response="pres")$mlik[1, 1])
    )

mod_comp_res <- bind_rows(
    mutate(mod_comp1, jsdm="independent residuals"), mutate(mod_comp2, jsdm="joint species residuals")
)

best_model <- slice_max(mod_comp_res, mlik, n=1)$form[[1]]

ggplot(mod_comp_res, aes(fx, rx)) +
    geom_tile(aes(fill=mlik), col="white") +
    geom_tile(fill=NA, col="orange", data=slice_max(mod_comp_res, mlik, n=1), size=1.3) +
    geom_text(aes(label=round(mlik, 1)), col="white", size=4.9) +
    facet_grid(.~jsdm, scales="free_x", space="free") +
    scale_fill_viridis_c(option="viridis", values=scales::rescale(sort(mod_comp_res$mlik))^3) +
    labs(x="Fixed effects", y="Random effects", fill="Marginal likelihood") +
    theme_bw() +
    scale_x_discrete(expand=expansion()) +
    scale_y_discrete(expand=expansion()) +
    theme(panel.spacing=unit(0,"lines"))

ggsave("figs/pres-model-performance.pdf", width=7.5, height=6)

# Best model performance----------------------------------------------------------

best_fit <- fit_model(best_model, parks_model_data)

## Marginal posteriors for fixed effect coefficients
best_fit$marginals.fixed |> 
    summ_fx() |> 
    mutate(species=str_extract(var, "Amblyomma americanum|Dermacentor variabilis|Ixodes spp\\."), var=fx_labels(var)) |> 
    ggplot(aes(m, var)) +
    geom_vline(xintercept = 0, col = 'gray55', linetype = 'dashed') +
    geom_linerange(aes(xmin = lo95, xmax = hi95), size = 1, col = 'darkorange2') +
    geom_linerange(aes(xmin = lo5, xmax = hi5), size = 2, col = 'darkorange2') +
    geom_point(size = 1.5, col = 'dodgerblue3') +
    facet_wrap(~species, nrow=3) +
    labs(x = 'Log odds', y = NULL) +
    theme_bw()

ggsave("figs/pres-fixed-effects.pdf", width=5, height=7)

## ROC curves for all observations and each species, based on mean fitted values
preds <- transmute(parks_model_data, tick_class, pres, pred=best_fit$summary.fitted.values$mean)

roc <- list(
    Dermacentor=roc(filter(preds, tick_class == "Dermacentor variabilis"), pres, pred),
    Ixodes=roc(filter(preds, tick_class == "Ixodes spp."), pres, pred),
    Amblyomma=roc(filter(preds, tick_class == "Amblyomma americanum"), pres, pred),
    All=roc(preds, pres, pred)
)

names(roc) <- imap_chr(roc, ~paste0(.y, " (AUC=", round(.x$auc, 4), ")"))

ggroc(roc, alpha=0.75, size=1.3) +
    theme_bw() +
    labs(x="Specificity", y="Sensitivity", col=NULL) +
    theme(
        legend.position=c(0.7, 0.2), 
        legend.background=element_rect(color="gray40"), 
        legend.text=element_text(size=rel(0.6))
    )

ggsave("figs/pres-roc.pdf", width=4, height=4)

## TODO: Temporal and spatial posterior trends
as_tibble(fit_pres$summary.random$`jan_min_temp`) |> 
    mutate(group=rep(1:3, each=102%/%3)) |> 
    ggplot(aes(ID, `0.5quant`)) +
    geom_line(aes(col=group, group=group)) +
    geom_point(data=transmute(pres, ID=jan_min_temp, `0.5quant`=pres))

# Make data for modeling abundance. All life_stages aggregated for this one
# Remove min_temp since corr. with max_temp and mean_rh and normalize
abun <- parks |> 
    filter(
        tick_class %in% c("Amblyomma americanum", "Ixodes spp."),
        life_stage != "larva"
    ) |> 
    group_by(date, site, tick_class) |>
    mutate(count=sum(count)) |> 
    ungroup() |> 
    select(-life_stage, -genus, -species) |> 
    distinct() |> 
    mutate(id=1:n())

terms <- abun |> 
    st_drop_geometry() |> 
    select(land_cover:mean_rh) |> 
    colnames() |> 
    str_c(collapse="+")

# Baseline regression-------------------------------------------------------------
fit_base <- inla(
    as.formula(paste0("count ~ 0 + tick_class + (", terms, "):tick_class")),
    data=st_drop_geometry(abun),
    family="nbinomial",
    control.fixed = list(
        expand.factor.strategy="inla",
        prec=0.25 #TODO: probs to small for counts data
    ),
    control.compute=list(dic=TRUE),
    control.predictor=list(link=1, compute=TRUE)
)

summary(fit_base)
compare_pred_counts(fit_base)
# Remember odd residuals come up here! So probably not just due to issues with fit; must be with covariates
plot_residuals(fit_base, abun) +
    ylim(-100, 100)

bad_res <- fit_base$summary.fitted.values |> 
    as_tibble() |>
    mutate(idx=1:n()) |> 
    slice_max(sd, n=8)

abun[bad_res$idx,]

# Abundance with B-spline---------------------------------------------------------

terms <- abun |> 
    st_drop_geometry() |> 
    select(tree_canopy:mean_rh) |> 
    colnames()

hyper_spline <- "list(theta=list(prior='pc.prec', param=c(1,0.01)))"
terms <- str_c("f(inla.group(", terms, "), model='rw1', hyper=", hyper_spline, ")", collapse="+")

fit_spline <- inla(
    as.formula(paste0("count ~ 0 + tick_class + land_cover +", terms)),
    # count ~ 0 + tick_class + land_cover + f(jan_min_temp, model = "rw1"),
    data=st_drop_geometry(abun),
    family="nbinomial",
    control.fixed = list(
        expand.factor.strategy="inla",
        prec=0.25
    ),
    control.compute=list(dic=TRUE),
    control.predictor=list(link=1, compute=TRUE)
)

summary(fit_spline)

as_tibble(fit_spline$summary.random$`inla.group(max_temp)`) |> 
    ggplot(aes(ID, mean)) +
    geom_line() +
    geom_point(data=transmute(abun, ID=max_temp, mean=log1p(count)))

##

fit_sep <- inla(
    count_bin ~ 0 + mean_temp:group + precipitation:group + mean_rh:group + 
        f(group, model="iid") + f(site, model="iid") + f(month, model="ar1"),
    data=parks,
    family="binomial",
    control.compute=list(return.marginals.predictor=TRUE, dic=TRUE)
)

summary(fit_sep)
compare_pred_counts(fit_sep)
plot_residuals(fit_sep)
plot_group_marginals(fit_sep)

# Maybe makes sense that group intercepts are bad here? I'm thinking they typically are not used, but check
# lit for any comments on this
fit_mon_jsdm <- inla(
    count_bin ~  0 + group + mean_temp:group + precipitation:group + mean_rh:group +
        # f(month, model="ar1") +
        # f(site, model="iid") +
        f(
            id, model="iid3d", n=nrow(parks),
            hyper=list(prec1=list(param=c(4, rep(1, 3), rep(0, 3))))
        ),
    data=parks, 
    family="binomial",
    control.compute=list(return.marginals.predictor=TRUE, dic=TRUE)
    # control.inla=list(int.strategy="ccd")
)

# When I changed to ZIP version 0, the correlations completely changed! That means zero counts were
# necessary to create the correlation pattern (implied by version 1)

# When I used neg binom, correlation changed/went away, residuals looked weird, I could decrease hyper, 
# and ML was better but DIC was worse (???)

inla.rerun(fit_mon_jsdm)

fit_mon_jsdm$summary.random$id

summary(fit_mon_jsdm)
compare_pred_counts(fit_mon_jsdm)
plot_residuals(fit_mon_jsdm)
plot_group_marginals(fit_mon_jsdm)

fit_mon_jsdm$summary.random$site

# TODO had a thought that some of this may be due to very non-gaussian correlations: ixodes and amblyomma tended
# to be completely exclusive, for example
as_tibble(fit_mon_giid$summary.random$month) |> 
    ggplot(aes(as.factor(ID), `0.5quant`)) +
    geom_point() +
    geom_errorbar(aes(ymin=`0.025quant`, ymax=`0.975quant`), width=0.2)

parks |> 
    mutate(count_pred=fit_mon_giid$summary.fitted.values$`0.5quant`) |> 
    ggplot(aes(date, log1p(count_pred), col=group)) +
    geom_point(aes(y=log1p(count)), col="black") +
    geom_line(shape=1, alpha=0.75)

plot_group_marginals(fit_mon_giid$marginals.random$group)
plot_group_marginals(fit_mon_giid$marginals.random$ggg)

qplot(
    filter(parks, group == "Amblyomma_larva")$count,
    filter(parks, group == "Ixodes_adult")$count
)

fit_mg_inter <- inla(
    count_bin ~  0 + mean_temp:group + precipitation:group + mean_rh:group + group +
            f(
            month, model="ar1",
            group=group_id, ngroup=3,
            control.group=list(model="ar1")
        ),
    #     f(
    #     month, model="rw1", cyclic=FALSE,
    #     group=group_id, ngroup=3,
    #     control.group=list(model="rw1", cyclic=FALSE)
    # ) +
    #     f(id, model="iid", hyper=list(prec=list(prior = "loggamma", param = c(1, 0.1)))),
    data=mutate(parks, group_id=as.numeric(fct_inorder(group))), 
    family="binomial",
    control.compute=list(return.marginals.predictor=TRUE, dic=TRUE)
    # control.inla=list(int.strategy="ccd", strategy="laplace")
)

inla.hyperpar(fit_mg_inter)

summary(fit_mg_inter)
compare_pred_counts(fit_mg_inter)
plot_residuals(fit_mg_inter)
plot_group_marginals(fit_mg_inter)

fit_mg_inter$summary.random$id

as_tibble(fit_mg_inter$summary.random$month) |>
    mutate(group=rep(unique(parks$group), each=8)) |> 
    ggplot(aes(as.factor(ID), `0.5quant`, col=group, group=1)) +
    geom_point(size=2) +
    geom_line() +
    geom_errorbar(aes(ymin=`0.025quant`, ymax=`0.975quant`), width=0.2, alpha=0.5) +
    facet_wrap(~group, scales="fixed") +
    labs(x="Site effect (log)", y="Month") +
    theme(legend.position="none")

# compare_pred_counts <- function(fit) {
#     parks |> 
#         mutate(count_pred=fit$summary.fitted.values$`0.5quant`) |> 
#         ggplot(aes(count, count_pred, col=group)) +
#         geom_point() +
#         geom_abline(slope=1, intercept=0, col="gray50")
# }
# 
# plot_residuals <- function(fit, data, x="id", col_by="tick_class") {
#     res <- data |> 
#         bind_cols(as_tibble(fit$summary.fitted.values)) |> 
#         mutate(rmean=count - `0.5quant`, rlow=count - `0.025quant`, rhigh=count - `0.975quant`)
#     
#     ggplot(res, aes(.data[[x]], col=.data[[col_by]])) +
#         geom_errorbar(aes(ymin=rlow, ymax=rhigh), width=0.4) +
#         geom_point(aes(y=rmean))
# }
# 
# plot_group_marginals <- function(fit, count_scale=FALSE) {
#     if (count_scale)
#         marg <- map(fit$marginals.random$group, ~inla.tmarginal(exp, .x))
#     else
#         marg <- fit$marginals.random$group
#     
#     group_marg <- map2_dfr(marg, unique(parks$group), ~tibble(x=.x[,1], y=.x[,2], group=.y))
#     
#     ggplot(group_marg, aes(x, y, col=group)) +
#         geom_line()
# }
