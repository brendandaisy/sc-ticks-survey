library(INLA)
library(tidyverse)
library(sf)

# Prepare a list of variables (just coefficients [inc. intercept] and not residuals right now) 
# for the "selection" option for INLA. Used to compute covariance matrix
sel_list_inla <- function(df) {
    vars <- c(
        str_c("land_cover", levels(df$land_cover)),
        "tree_canopy", "elevation", "jan_min_temp", "max_temp", "precipitation", "mean_rh"
    )
    coefs <- flatten_chr(map(
        vars,
        ~str_c("tick_class", c("Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp."), ":", .x)
    ))
    coefs <- c(str_c("tick_class", c("Amblyomma americanum", "Dermacentor variabilis", "Ixodes spp.")), coefs)
    sel <- rep(1, length(coefs)) |> as.list()
    names(sel) <- coefs
    return(sel)
}

# returns a matrix of risk probabilities with locations as rows and `n` cols
rpost_predict <- function(obs_df, pred_idxs, n=1) {
    f <- formula_jsdm(obs_df)
    ft <- fit_model(f, obs_df, fx_prec=0.2, selection=list(Predictor=pred_idxs))
    jsamp <- inla.rjmarginal(n, ft$selection)
    return(inla.link.invlogit(jsamp$samples))
}

# Bayesian D-optimality criteria
# assumes new_df has sampled pres, with parks data added
util_bd_rep <- function(new_df, sel, stable=TRUE) {
    f <- formula_jsdm(new_df)
    ft <- fit_model(f, new_df, fx_prec=0.2, selection=sel)
    -log(det(ft$selection$cov.matrix))
}

# optional full_df can be used to extract all info if d_df only has site and date
# this can be useful since there are J rows per visit
utility <- function(d_df, known_df, n=1, full_df=NULL, by=c("date", "site"), u_only=FALSE) {
    if (!is.null(full_df)) {
        d_df <- semi_join(full_df, d_df, by=by)
    }
    # setup new data frame (parks + d_df) and get prediction matrix
    sel <- sel_list_inla(known_df)
    new_df <- prep_new_data(known_df, d_df, scale=FALSE)
    pred_idxs <- which(is.na(new_df$pres))
    pred_risk_mat <- rpost_predict(new_df, pred_idxs, n=n)
    # for each row in pred mat, get a sample from p(y_pred|y) and calc util
    u_rep <- map_dbl(1:ncol(pred_risk_mat), ~{
        new_df$pres[pred_idxs] <- rbinom(nrow(d_df), rep(1, nrow(d_df)), pred_risk_mat[,.x])
        util_bd_rep(new_df, sel)
    })
    if (u_only)
        return(mean(u_rep))
    return(tibble_row(utility=mean(u_rep), design=list(d_df)))
}

# if file exists, it will be appended
# the result is returned
save_util_res <- function(
        utils, strat=c("random", "simple-var", "bayes-opt", "sim-ann"), n, 
        append=TRUE, alpha=NULL
) {
    a <- if (is.null(alpha)) "" else paste0("-alpha=", alpha)
    file <- paste0("util-exper/", strat, "-n=", n, a, ".rds")
    if (file.exists(file))
        res <- readRDS(file)
    else
        res <- tibble()
    ret <- bind_rows(res, utils) |> mutate(strat=strat)
    saveRDS(ret, file)
    return(ret)
}