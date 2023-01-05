library(INLA)
library(tidyverse)
library(sf)
# library(MASS, exclude=c("select"))
# library(furrr)
# library(gridExtra)
# library(DiceOptim)
library(fields)

source("other-helpers.R")
source("utility-helpers.R")

# Add columns for data points in normalized design space, $(X, Y, t) \in [0, 1]^3$
# Used for Bayesian Optimization
sp_add_norm <- function(df) {
    ret <- df |> 
        st_normalize()
    
    ret |> 
        st_drop_geometry() |> 
        mutate(X=st_coordinates(ret)[,1], Y=st_coordinates(ret)[,2], t=(month-3) / 9, .before=tick_class)
}

# Data preparation----------------------------------------------------------------
parks_model_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

all_data <- append_pred_data(parks_model_data) |> 
    sp_add_norm()

# Predict expected risk across grid, for sampling---------------------------------

# fit_all <- fit_model(formula_jsdm(all_data), all_data, fx_prec=0.2)
fit_all <- readRDS("inla-fit-all-data.rds")
fit_all$summary.fixed # assert: fixed effects were included in the model correctly?

post_pred_all <- all_data |>
    mutate(
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd,
        var_eta=fit_all$summary.linear.predictor$sd^2
    )

pp_grid <- filter(post_pred_all, is.na(pres))
pp_parks <- filter(post_pred_all, !is.na(pres))

# Utility functions considered----------------------------------------------------

B <- seq(5, 30, 5)
reps <- 5
expg <- expand_grid(B=B, rep=1:reps)
n <- 100

# Random selection, baseline design strategy--------------------------------------
u_random <- pmap_dfr(expg, ~{
    d <- pp_grid |> 
        select(date, site) |> 
        slice_sample(n=..1)
    
    mutate(util_bd_opt(d, pp_parks, n=n, full_df=pp_grid), B=..1)
})

u_random <- save_util_res(u_random, "random", n=n)

# Select sites with largest average variance among species,-----------------------
# only choosing a site one time---------------------------------------------------
pp_grid_sort <- pp_grid |> 
    group_by(date, site) |> 
    # average variance from each species:
    summarise(var_eta=mean(var_eta), .groups="drop") |> 
    arrange(desc(var_eta)) |> 
    group_by(site) |> 
    # only allow each site to be selected:
    filter(row_number() == 1) |> 
    ungroup()

u_simple_var <- imap_dfr(unique(expg$B), ~{ # only have to do each B once since deterministic
    print(paste0("i=", .y))
    d <- head(pp_grid_sort, .x)
    mutate(util_bd_opt(d, pp_parks, n=n, full_df=pp_grid), B=.x)
})

u_simple_var <- save_util_res(u_simple_var, "simple-var", n=n)

# TODO: makes sense! Next is to do regular spacing, or choosing locations as far away as possible perhaps
# also consider something slightly more exploratory but with simple variance as based
# e.g. choose 100 new locations and pick one with top variance

d_resamp <- filter(pp_parks, lubridate::year(date) == 2021)
u_resamp <- util_bd_opt(d_resamp, pp_parks, n=n)
u_none <- util_bd_opt_rep(prep_new_data(pp_parks, scale=FALSE), sel_list_inla(pp_parks))

slice_sample(pp_grid, n=10) |> 
    mutate(land_cover=as.numeric(land_cover)) |> 
    select(land_cover:mean_rh) |>
    dist()

##TODO do this for the normalize (x, y, t) space, since already should have code for this

geo_grid <- all_data |> 
    filter(str_detect(site, "g\\d+")) |> 
    select(site) |> 
    distinct()

geo_parks <- all_data |> 
    filter(!str_detect(site, "g\\d+")) |> 
    select(site) |> 
    distinct()

# Bayesian optimization-----------------------------------------------------------

init_utils <- function(pred_df, known_df, n_init=10, n=10) {
    eval_pts <- pred_df |> 
        select(X, Y, t) |> 
        distinct() |> 
        slice_sample(n=n_init)
    
    map_dfr(1:n_init, ~{
        d <- eval_pts[.x,]
        tibble_row(
            X=d$X, Y=d$Y, t=d$t, 
            utility=util_bd_opt(d, known_df, n=n, full_df=pred_df, tibble=FALSE, by=c("X", "Y", "t"))
        )
    })
}

iter_bayes_opt <- function(
        init_utils, pred_df, known_df, d_so_far=tibble(), n=10,
        max_iter=10, max_same_iter=max_iter, noise_util=10, pick_top=5,
        control_km=list(trace=FALSE, pop.size=100, max.generations=50, wait.generations=10),
        control_akg=list(pop.size=100, max.generations=50, wait.generations=10, BFGSburnin=10)
) {
    init_utils$utility <- -init_utils$utility # make sure is minimization problem
    best_f <- min(init_utils$utility)
    nsame <- 0
    # setup record of considered designs, with blank info
    ures <- bind_rows(init_utils, tibble(utility=rep(NA, max_iter)))
    
    # eval_pts <- pred_loc |>
    #     distinct(X, Y, t) |> 
    #     relocate(X, Y, t)
    
    for (i in seq_len(max_iter)) { # choose a new point to evaluate U
        ucur <- ures |> filter(!is.na(utility)) # current design surface
        krig <- km(
            ~1, select(ucur, X, Y, t), ucur$utility, scaling=TRUE,
            noise.var=rep(noise_util, nrow(ucur)), lower=rep(0.1, 3), upper=rep(1, 3),
            optim.method="BFGS"
        )
        return(krig)
        res <- max_AKG(
            krig, new.noise.var=noise_util, type="UK", lower=rep(0, 3), upper=rep(1, 3),
            parinit=unlist(select(slice_sample(ucur, n=1), -utility)),
            control=control_akg
        )
        print(res)
        next_pt <- tibble_row(X=res$par[1], Y=res$par[2], t=res$par[3])
        d_cand <- bind_rows(d_so_far, next_pt)
        next_u_eval <- util_bd_opt(d_cand, known_df, n=n, full_df=pred_df, tibble=FALSE, by=c("X", "Y", "t"))
        ures[nrow(init_utils)+i,] <- mutate(next_pt, utility=-next_u_eval)
        # sub_eval_pts <- slice_sample(eval_pts, n=2500)
        # eqi_pts <- map_dbl(
        #     1:nrow(sub_eval_pts), 
        #     ~EQI(c(sub_eval_pts[.x,]), krig, new.noise.var=noise_util, type="UK")
        # )
        # next_pt <- sub_eval_pts[which.max(eqi_pts),]
        # next_u_eval <- util_bd_opt(next_pt, )

        if (res$value <= best_f) {
            best_f <- res$value
            nsame <- 0
        } else {
            nsame <- nsame + 1
        }
        if (nsame > max_same_iter)
            break
    }
    ures <- filter(ures, !is.na(utility))
    
    return(list(
        best=slice_min(ures, utility, n=pick_top, with_ties=FALSE),
        all=ures
    ))
}

utils0 <- init_utils(pp_grid, pp_parks, n_init=30, n=5)

bo1 <- iter_bayes_opt(utils0, pp_grid, pp_parks, n=5)

# Simulated Annealing-------------------------------------------------------------

# T0 is initial order of mag. worse a prop can be, for acc. prob. 0.37
# alpha > 1 means faster cooling (less acceptance of worse moves)
simulated_annealing <- function(
        pred_df, known_df, B=5, iter=10, alpha=1, T0=1, n_est=5, 
        d_so_far=tibble(), d_init=NULL, weight_by=var_eta
) {
    d_cand <- pred_df |> 
        group_by(date, site) |> 
        summarise(var_eta=mean(var_eta), .groups="drop")
    
    if (is.null(d_init))
        d_curr <- d_best <- slice_sample(d_cand, n=B, weight_by={{weight_by}}) |> bind_rows(d_so_far)
    else
        d_curr <- d_best <- bind_rows(d_so_far, d_init)
    u_curr <- u_best <- util_bd_opt(d_curr, known_df, n=n_est, full_df=pp_grid, tibble=FALSE)
    i_best <- 0
    
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    for (i in 1:iter) {
        next_pt <- anti_join(d_cand, d_curr, by=c("date", "site")) |> 
            slice_sample(n=1, weight_by={{weight_by}})
        
        d_prop <- slice_sample(d_curr, n=nrow(d_curr)-1) |> 
            bind_rows(next_pt)
        
        u_prop <- util_bd_opt(d_prop, pp_parks, n=n_est, full_df=pp_grid, tibble=FALSE)
        if (u_prop > u_best) {
            i_best <- i
            d_curr <- d_best <- d_prop
            u_curr <- u_best <- u_prop
        } else if (u_prop > u_curr) {
            d_curr <- d_prop
            u_curr <- u_prop
        } else if (runif(1) < exp((log10(u_prop) - log10(u_curr))/T_sched[i])) {
            d_curr <- d_prop
            u_curr <- u_prop
        }
        print(paste0("Current log(U)=", round(log(u_curr), 3), " at i=", i))
    }
    return(list(d=d_best, u_approx=u_best, at=i_best))
}

sa_res <- simulated_annealing(pp_grid, pp_parks, B=5, iter=120, alpha=0.9, n_est=50, d_so_far=sa_res$d)

u_sim_ann <- util_bd_opt(sa_res$d, pp_parks, n=n, full_df=pp_grid) |> 
    mutate(B=nrow(sa_res$d))

saveRDS(sa_res, "util-exper/sa-test-alpha=0.9-nest=15-iter=200.rds")

save_util_res(u_sim_ann, "sim-ann", n=n)

# TODO clear this works now. Plan I'd say is to write methods now.
# After any minor updates that become known I want to do, its time to
# run the non-random strategies with like n=200 and LEAVE THINGS AS FINAL

# Plot the results comparing all strategies---------------------------------------

bind_rows(
    mutate(readRDS("util-exper/random-n=50.rds"), strat="random"),
    mutate(u_simple_var, strat="simple-var"),
    u_sim_ann
) |> 
    ggplot(aes(B, utility, col=strat)) +
    geom_point() +
    # geom_hline(yintercept=u_resamp$utility, col="green") +
    # geom_hline(yintercept=u_none, col="gray40")
    scale_y_log10() +
    annotation_logticks(sides="l")
