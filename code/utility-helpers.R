# --------------------------------------------------------------------------------
# utility-helpers.R---------------------------------------------------------------
# functions for calculating utility and for the two design criteria used----------
# --------------------------------------------------------------------------------
library(INLA)
library(tidyverse)
library(sf)
library(furrr)

#' @name utility
#' @description Workhorse function to compute the utility of a proposed design `d`
#'
#' @param d A df of visits, containing at least columns `date` and `site`
#' @param known_df A df of any initial collections data
#' @param n The number of Monte Carlo samples used to approximate utility
#' @param pred_df Optional df containing info for proposed visits, to retrieve e.g.
#' covariates if `d` only contains two columns. Probably output from `append_pred_data`
#' @param util_fun Function defining design criterion
#' @param u_only Whether to produve a tibble containing results or just a number
#' @param copy For if `inner_join` is giving a "different data source" error
#' @param by Passed to `inner_join` to silence warning
#' @param sel Passed to `util_fun`, used for `util_dopt`
#' @param risk_df Passed to `util_fun`, used for `util_risk_sd`
#'
#' @return A tibble or double of the utility
utility <- function(
        d, known_df, n=1, pred_df=NULL, util_fun=util_dopt, u_only=FALSE,
        copy=FALSE, by=c("date", "site"), sel=NULL, risk_df=NULL
) {
    # grab the covariates etc. from `pred_df` if necessary
    d_df <- if (!is.null(pred_df)) inner_join(pred_df, d, by=by, copy=copy) else d
    # combine initial and proposed visits and do posterior prediction
    new_df <- bind_rows(known_df, d_df)
    pred_idxs <- which(is.na(new_df$pres))
    pred_risk_mat <- rpost_predict(new_df, pred_idxs, n=n)
    
    # for each row in pred mat, get a sample from P(y_pred|y) and calc utility
    u_rep <- future_map_dbl(1:ncol(pred_risk_mat), ~{
        new_df$pres[pred_idxs] <- rbinom(nrow(d_df), rep(1, nrow(d_df)), pred_risk_mat[,.x])
        tryCatch(
            util_fun(new_df, sel=sel, risk_df=risk_df),
            error=function(e) {on_inla_error(e, d, new_df$pres[pred_idxs])}
        )
    }, .options=furrr_options(seed=TRUE)
    )
    # format and return the results, averaging U(d, y) over y
    if (u_only)
        return(mean(u_rep, na.rm=TRUE))
    uname <- substitute(util_fun) |> as.character() |> str_remove("util_")
    return(tibble_row({{uname}} := mean(u_rep, na.rm=TRUE), design=list(d)))
}

#' @name rpost_predict
#' Fit an INLA model to sample presence/absence for proposed future visits 
#'
#' @param new_df Initial visitation data
#' @param pred_idxs The indices/rows in `new_df` corresp. to future visits. Posterior prediction is done at these sites
#' @param n The number of Monte Carlo samples used to approximate utility
#'
#' @return a matrix of risk probabilities with locations as rows and `n` cols
rpost_predict <- function(new_df, pred_idxs, n=1) {
    f <- formula_bed()
    ft <- fit_model(f, new_df, fx_prec=0.2, selection=list(Predictor=pred_idxs))
    jsamp <- inla.rjmarginal(n, ft$selection)
    return(inla.link.invlogit(jsamp$samples))
}

# Bayesian D-optimality criteria
# assumes new_df has sampled pres, with parks data added
util_dopt <- function(new_df, sel, risk_df=NULL) {
    f <- formula_bed()
    ft <- fit_model(f, new_df, fx_prec=0.2, selection=sel)
    -log(det(ft$selection$cov.matrix))
}

# Criteria 2 for maximizing info for high-risk visits
util_risk_sd <- function(new_df, risk_df, sel=NULL) {
    ndf <- bind_rows(new_df, risk_df)
    f <- formula_bed()
    ft <- fit_model(f, ndf, fx_prec=0.2)
    vardiff <- risk_df$sd_pres - ft$summary.fitted.values$sd[-c(1:nrow(new_df))]
    return(max(vardiff))
}

#' sel_list_inla
#' Prepare a list of the  environmental effects for the `selection` option for INLA
#' Used to compute the posterior covariance matrix used for D-opt criterion
#'
#' @param df To get the environmental variable names
#'
#' @return A named list
sel_list_inla <- function(df) {
    vars <- c(
        str_c("land_cover", levels(df$land_cover)),
        "tree_canopy", "elevation", "jan_min_temp", "max_temp", "precipitation", "mean_rh"
    )
    # add the tick-specific intercepts
    coefs <- c(str_c("tick_class", c("Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp.")), vars)
    sel <- rep(1, length(coefs)) |> as.list()
    names(sel) <- coefs
    return(sel)
}

on_inla_error <- function(e, design, sample) {
    print(e)
    print("Design:")
    print(as.data.frame(design))
    print("Sample:")
    print(sample)
    stop()
}
