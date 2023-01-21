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

inla.setOption(inla.mode="classic")

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
parks_data <- read_parks_sf("geo-files/parks-with-covars.shp", drop=min_temp) |> 
    prep_parks_model_data(rescale=FALSE) # don't rescale, since scaled wrt. grid data below

# Predict expected risk across grid, for sampling---------------------------------

# fit_all <- fit_model(formula_jsdm(all_data), all_data, fx_prec=0.2)
fit_all <- readRDS("inla-fit-all-data.rds")
fit_all$summary.fixed # assert: fixed effects were included in the model correctly?

all_data <- append_pred_grid(parks_data) |> 
    sp_add_norm() |> 
    mutate(
        mean_pres=fit_all$summary.fitted.values$mean,
        sd_pres=fit_all$summary.fitted.values$sd,
        var_eta=fit_all$summary.linear.predictor$sd^2
    )

grid_mod <- filter(all_data, is.na(pres)) # final not yet observed data for fitting models
parks_mod <- filter(all_data, !is.na(pres)) # final observed data for fitting models

# Utility functions considered----------------------------------------------------

num_loc <- c(1, seq(5, 20, 5))
reps <- 5
expg <- expand_grid(num_loc=num_loc, rep=1:reps)
n <- 20

# Random selection, baseline design strategy--------------------------------------
u_random <- pmap_dfr(expg, ~{
    d <- grid_mod |> 
        select(date, site) |> 
        slice_sample(n=..1)
    
    mutate(utility(d, parks_mod, n=n, full_df=grid_mod), num_loc=..1, rep=..2)
})

# u_random <- save_util_res(u_random, "random", n=n, append=FALSE)
u_random <- load_util_res(u_random, "random", n=n)

# Select sites with largest average variance among species,-----------------------
# only choosing a site one time---------------------------------------------------
grid_mod_sort <- grid_mod |> 
    group_by(date, site) |> 
    # average variance from each species:
    summarise(var_eta=mean(var_eta), .groups="drop") |> 
    arrange(desc(var_eta)) |> 
    group_by(site) |> 
    # only allow each site to be selected once:
    filter(row_number() == 1) |> 
    ungroup()

u_simple_var <- imap_dfr(num_loc, ~{ # only have to do each B once since deterministic
    print(paste0("i=", .y))
    d <- head(grid_mod_sort, .x)
    mutate(utility(d, parks_mod, n=n, full_df=grid_mod), num_loc=.x)
})

u_simple_var <- save_util_res(u_simple_var, "simple-var", n=n)
u_simple_var <- load_util_res(u_simple_var, "simple-var", n=n)

###

tmp <- grid_mod |> 
    select(date, site) |> 
    slice_sample(n=1)

tmp <- semi_join(grid_mod, tmp, by=c("date", "site")) |> 
    bind_rows(parks_mod)

X <- dummyVars(
    pres ~ -1 + tick_class + tick_class:land_cover + tick_class:tree_canopy + tick_class:elevation + 
        tick_class:jan_min_temp + tick_class:max_temp + tick_class:precipitation + tick_class:mean_rh, 
    data=mutate(tmp, land_cover=fct_drop(land_cover)),
    fullRank=TRUE
)
X <- predict(X, tmp)

get_eta <- function(x) {
    eta <- sum(x * fit_all$summary.fixed$mean)
    exp(eta) / (1 + exp(eta))^2
}

tmp2 <- map(1:nrow(X), ~X[.x,] %*% get_eta(Xfull[.x,]) %*% t(X[.x,]))

I <- t(X) %*% diag(tmp$sd_pres^2) %*% X
I <- reduce(tmp2, ~.x + .y)

eigen(I)$values < -1e8
cov_mat <- solve(I)

image(I)

I %*% cov_mat
det(cov_mat)
sum(diag(cov_mat))

# "One-off" design strategies-----------------------------------------------------

## Utility with no additional sampling
u_none <- util_rep(prep_new_data(parks_mod, scale=FALSE), sel_list_inla(parks_mod))

## Reuse the same schedule from 2021
d_resamp <- filter(parks_mod, lubridate::year(date) == 2021)
d_resamp$pres <- NA
u_resamp <- utility(d_resamp, parks_mod, n=n)

## Visit all 30 parks again in December
all_data_sp <- append_pred_grid(parks_data)

d_parks_dec <- parks_data |> 
    distinct(site, geometry) |> 
    st_join(filter(all_data_sp, month == 12), st_nearest_feature) |> 
    st_drop_geometry() |> 
    rename(site=site.x) |> 
    select(-site.y)

d_parks_dec$pres <- NA
u_parks_dec <- utility(d_parks_dec, parks_mod, n=n)

rm(all_data_sp) # remove this for space since no longer needed

# slice_sample(pp_grid, n=10) |> 
#     mutate(land_cover=as.numeric(land_cover)) |> 
#     select(land_cover:mean_rh) |>
#     dist()

# Bayesian optimization-----------------------------------------------------------

init_utils <- function(pred_df, known_df, d_so_far=tibble(), n_init=10, n=10) {
    eval_pts <- pred_df |> 
        distinct(X, Y, t) |> 
        anti_join(d_so_far) |> 
        slice_sample(n=n_init)
    
    map_dfr(1:n_init, ~{
        d <- eval_pts[.x,]
        print(d)
        tibble_row(
            X=d$X, Y=d$Y, t=d$t, 
            utility=utility(bind_rows(d_so_far, d), known_df, n=n, full_df=pred_df, u_only=TRUE, by=c("X", "Y", "t"))
            # full_df is everything even when d_so_far is not empty
        )
    })
}

iter_bayes_opt <- function(
        init_utils, pred_df, known_df, d_so_far=tibble(), n=10,
        max_iter=10, max_same_iter=max_iter, noise_util=0.1,
        control_km=list(trace=FALSE, pop.size=100, max.generations=50, wait.generations=10),
        control_akg=list(pop.size=100, max.generations=50, wait.generations=10, BFGSburnin=10)
) {
    init_utils$utility <- -init_utils$utility # make sure is minimization problem
    best_f <- min(init_utils$utility)
    nsame <- 0
    # setup record of considered designs, with blank info
    ures <- bind_rows(init_utils, tibble(utility=rep(NA, max_iter)))
    
    eval_pts <- pred_df |> # note that it is allowed to select a pt already in d_so_far (=visit 2x same month)
        distinct(X, Y, t) |>
        relocate(X, Y, t) # ensure they are in the right order for kriging
    
    for (i in seq_len(max_iter)) { # choose a new point to evaluate U
        ucur <- ures |> filter(!is.na(utility)) # current design surface
        
        krig <- km(
            ~1, select(ucur, X, Y, t), ucur$utility, scaling=TRUE,
            noise.var=rep(noise_util, nrow(ucur)), lower=rep(0.1, 3), upper=rep(1, 3),
            optim.method="BFGS", control=control_km
        )
        # res <- max_AKG(
        #     krig, new.noise.var=noise_util, type="UK", lower=rep(0, 3), upper=rep(1, 3),
        #     parinit=unlist(select(slice_sample(ucur, n=1), -utility)),
        #     control=control_akg
        # )
        # print(res)
        # next_pt <- tibble_row(X=res$par[1], Y=res$par[2], t=res$par[3])
        # d_cand <- bind_rows(d_so_far, next_pt)
        # next_u_eval <- util_bd_opt(d_cand, known_df, n=n, full_df=pred_df, tibble=FALSE, by=c("X", "Y", "t"))
        
        sub_eval_pts <- slice_sample(eval_pts, n=3000)
        eqi_pts <- map_dbl(
            1:nrow(sub_eval_pts),
            ~AKG(c(sub_eval_pts[.x,]), krig, new.noise.var=noise_util, type="SK")
        )
        next_pt <- sub_eval_pts[which.max(eqi_pts),]
        print(next_pt)
        d_cand <- bind_rows(d_so_far, next_pt)
        next_u_eval <- -utility(d_cand, known_df, n=n, full_df=pred_df, u_only=TRUE, by=c("X", "Y", "t"))
        print(next_u_eval)
        ures[nrow(init_utils)+i,] <- mutate(next_pt, utility=next_u_eval)
        print(as.data.frame(ures))

        if (next_u_eval <= best_f) {
            best_f <- next_u_eval
            nsame <- 0
        } else {
            nsame <- nsame + 1
        }
        if (nsame > max_same_iter)
            break
    }
    ures <- filter(ures, !is.na(utility)) # remove any unused rows from terminating early
    ures$utility <- -ures$utility # convert back to a maximization problem
    
    return(list(
        best=slice_max(ures, utility, n=1, with_ties=FALSE),
        all=ures
    ))
}

library(DiceKriging)
library(DiceOptim)

# TODO: these can be reused from d_random (but just for first step when num_loc=1)
utils1 <- init_utils(grid_mod, parks_mod, n_init=40, n=40)
bo1 <- iter_bayes_opt(utils1, grid_mod, parks_mod, n=40, max_iter=15)

utils2 <- init_utils(grid_mod, parks_mod, bo1$best, n_init=40, n=40)
bo2 <- iter_bayes_opt(utils2, grid_mod, parks_mod, d_so_far=bo1$best[1, c("X", "Y", "t")], n=40, max_iter=15)

utils3 <- init_utils(grid_mod, parks_mod, bind_rows(bo1$best, bo2$best), n_init=40, n=40)
bo3 <- iter_bayes_opt(bo3$all, grid_mod, parks_mod, d_so_far=bind_rows(bo1$best, bo2$best), n=40, max_iter=15)

ggplot(bo3$all, aes(X, Y, col=utility, size=utility)) +
    geom_point() +
    facet_wrap(~t, nrow=2)

utils0 |> 
    arrange(desc(utility)) |> 
    select(X:t) |> 
    slice(1) |> 
    utility(parks_mod, n=1, full_df=grid_mod, by=c("X", "Y", "t"))

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

u_one_offs <- tibble(
    strat=c("None (0)", "Repeat 2021 schedule (111)", "Revisit in December (30)"), 
    utility=c(u_none[1], u_resamp$Dfixed, u_parks_dec$Dfixed)
)

bind_rows(u_random, u_simple_var) |>
    ggplot(aes(num_loc, utility, col=strat)) +
    # geom_hline(aes(yintercept=utility, col=strat), data=u_one_offs, show.legend=FALSE) +
    geom_hline(yintercept=u_one_offs$utility, col=c("gray40", "blue", "green")) +
    annotate("text", label=u_one_offs$strat, y=u_one_offs$utility, col=c("gray40", "blue", "green"), x=22.5, hjust="left") +
    geom_point() +
    # geom_hline(yintercept=bo1$best$utility[1], col="red") +
    coord_cartesian(xlim = c(0, 21), clip = "off") +
    labs(x="Number of new visits", y="Utility", col="Search strategy") +
    theme_bw() +
    theme(plot.margin = unit(c(0, 1, 0, 0), "in"))

ggsave("figs/utility-results.pdf", width=5.3, height=3.8)
