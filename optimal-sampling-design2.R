library(INLA)
library(tidyverse)
library(sf)
# library(MASS, exclude=c("select"))
# library(furrr)
# library(gridExtra)
# library(DiceOptim)
library(fields)
library(caret)

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

# fit_all <- fit_model(formula_bed(), all_data, fx_prec=0.2)
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

u_simple_var <- map_dfr(num_loc, ~{ # only have to do each B once since deterministic
    print(paste0("num_loc=", .x))
    d <- head(grid_mod_sort, .x)
    mutate(utility(d, parks_mod, n=50, full_df=grid_mod), num_loc=.x)
})

saveRDS(u_simple_var, "util-results/util-simple-var.rds")

u_simple_var <- save_util_res(u_simple_var, "simple-var", n=n)
u_simple_var <- load_util_res(u_simple_var, "simple-var", n=n)

###

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

## 1. Utility with no additional sampling
u_none <- util_rep(obs_mod, sel_list_inla(obs_mod))

## 2-3. Visit all 30 parks again in December/June

d_parks_dec <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-12-01"))
    
u_parks_dec <- utility(d_parks_dec, obs_mod, n=n, pred_df=pred_mod)

d_parks_jun <- obs_mod |> 
    distinct(site) |> 
    mutate(date=as.Date("2023-06-01"))

u_parks_jun <- utility(d_parks_jun, obs_mod, n=n, pred_df=pred_mod)

## 4. Reuse the same schedule from 2021
d_resamp <- filter(obs_mod, lubridate::year(date) == 2021)
d_resamp$pres <- NA
u_resamp <- utility(d_resamp, obs_mod, n=n)

u_one_offs <- tibble(
    strat=c("    None (0)", "    Revisit in\n    June (30)", "    Revisit in\n    December (30)", "    Repeat 2021\n    schedule (111)"), 
    design=list(tibble(), distinct(d_parks_jun, date, site), distinct(d_parks_dec, date, site), distinct(d_resamp, date, site)),
    utility=c(u_none[1], u_parks_jun$dopt, u_parks_dec$dopt, u_resamp$dopt)
)

saveRDS(u_one_offs, "util-results/util-one-offs.rds")

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


# Plot the results comparing all strategies---------------------------------------

strats <- c("util-random", "util-simple-var", "util-sim-ann-local-alpha0pt9-T02", "util-sim-ann")
names(strats) <- c("Random", "Simple Var", "Local", "Simulated Annealing")

u_res <- strats |> 
    map_dfr(~mutate(readRDS(paste0("util-results/", .x, ".rds")), strat=.x)) |> 
    mutate(strat=fct_recode(strat, !!!strats))

u_res_random <- filter(u_res, strat == "Random")

u_res <- u_res_random |> 
    filter(num_loc == 1) |> 
    slice_max(dopt, n=1) |> 
    mutate(strat="Simulated Annealing") |> 
    bind_rows(u_res)

u_rand_mean <- u_res_random |> 
    group_by(num_loc) |> 
    summarise(mean=mean(dopt))
    
gg <- ggplot(res, aes(factor(num_loc), dopt)) +
    geom_violin(fill="#9ac9e7", col="gray70", alpha=0.75)  +
    geom_point(col="#9ac9e7", size=0.9) +
    coord_cartesian(clip="off") +
    # ylim(NA, 30) +
    labs(x="Number of new visits", y="Utility", col="Search strategy") +
    theme_bw() +
    theme(
        plot.margin=unit(c(t=0.3, r=1, b=0, l=0), "in"), 
        # axis.title=element_text(size=rel(1.1)),
        # axis.text=element_text(size=rel(1.1))
    )

ggsave("figs/util-results1.pdf", width=4.8, height=4)

hline_col <- c("#3c4276", "#3798a9", "#7fd7b8", "#98d180")

plot_oneoffs <- function(i) {
    gg2 <- gg +
        geom_hline(yintercept=u_one_offs$utility[1:i], col=hline_col[1:i], alpha=0.65, linewidth=1, linetype="dashed") +
        annotate(
            "text", 
            label=u_one_offs$strat[1:i], 
            y=u_one_offs$utility[1:i] + c(-0.2, -0.63, 0.63, 0)[1:i], 
            col=hline_col[1:i], x=5.5, hjust="left", size=3
        )
    
    ggsave(paste0("figs/util-results2", letters[i], ".pdf"), width=4.8, height=4)
    gg2
}

gg2 <- plot_oneoffs(4)

plot_searches <- function(vars) {
    u_sub <- filter(u_res, strat != "Random", strat %in% vars)
    
    gg2 +
        geom_line(aes(col=fct_relevel(strat, "Simple Var"), group=strat), u_sub) +
        geom_point(aes(col=fct_relevel(strat, "Simple Var")), u_sub, size=0.9) +
        scale_color_manual(values=c("#f5b43d", "#f53db5", "#c63df5")[seq_along(vars)]) +
        theme(legend.position="top", plot.margin=unit(c(t=0, r=1, b=0, l=0), "in"))
    
    ggsave(paste0("figs/util-results3", letters[length(vars)], ".pdf"), width=4.9, height=4.2)
}

plot_searches(c("Simple Var", "Local", "Simulated Annealing"))

###

all_loc_sp <- append_pred_grid(parks_data) |> 
    distinct(site, geometry)

plot_design_map <- function(d, col) {
    d_sp <- inner_join(all_loc_sp, d) |> 
        mutate(month=fct_drop(lubridate::month(date, label=TRUE, abbr=FALSE), c("January", "February")))
    
    ggplot(d_sp) +
        geom_sf(data=sc_state, fill=NA, col="gray70", linewidth=0.7, alpha=0.9) +
        geom_sf(col=col, alpha=0.8, size=1.2) +
        facet_wrap(~month, drop=FALSE, nrow=3) +
        map_theme
}

sc_state <- st_read("geo-files/south-carolina-county-boundaries.shp") |> 
    st_transform(st_crs(all_loc_sp)) |> 
    st_union() |> 
    nngeo::st_remove_holes()

map_theme <- theme_bw() +
    theme(
        axis.text=element_blank(), axis.ticks=element_blank(),
        panel.spacing=unit(0.55, "mm"),
        plot.margin=unit(c(0, 0.4, 0, 0), "mm"),
        panel.grid=element_blank(),
        strip.text=element_text(size=rel(1.15))
    )

plot_design_map(u_one_offs$design[[3]], hline_col[3])
ggsave("figs/util-map2c.pdf", width=3.8, height=3.2)

search_cols <- c("#f5b43d", "#f53db5", "#c63df5")
plot_design_map(filter(u_res, num_loc == 20, strat == "Simulated Annealing")$design[[1]], search_cols[3])
ggsave("figs/util-map3c.pdf", width=3.8, height=3.2)
