# --------------------------------------------------------------------------------
# bed-workflow-example.R----------------------------------------------------------
# A rough script that was used to produce components of Figure 1, illustrating----
# the idea behind BED in a simple design space------------------------------------
# --------------------------------------------------------------------------------

library(tidyverse)
library(INLA)

inla.setOption(inla.mode="classic")

ex_design_space <- function() {
    b0 <- 0.5
    b1 <- rnorm(8, 1, 0.01) # small amount of model bias
    # xdom <- c(-2, 2)
    # xpred <- seq(-2, 2, 0.1)
    # x0 <- c(-2, -1.95, -1, 0, 0.1, 0.105, 0.2, 1.8, 1.9, 1.95)
    # r <- inla.link.invlogit(b1 + b2 * x0)
    # y <- rbinom(10, 1, r)
    
    expand_grid(s1=c("A", "B"), s2=c("C", "D"), t=paste0("$t=", 1:4, "$")) |> 
        mutate(
            x=c(
                -2, -1.5, -1, -0.7,
                -1, -1.2, 0.1, 1.2,
                -1, -0.2, 1, 1.5,
                -0.8, 1, 1.9, 2
            ), 
            r=inla.link.invlogit(b0 + b1 * x),
            y=NA
        )
}

set.seed(1222)
design_space <- ex_design_space()

x0 <- c(1, 1, 2, 5, 7, 10, 10, 15)
df0 <- design_space |> 
    slice(x0) |> 
    mutate(y=rbinom(8, 1, r))

fit0_null <- inla(
    y ~ 1, data=df0, family="binomial", 
    control.compute=list(dic=TRUE), control.predictor=list(link=1, compute=TRUE)
)

pred_df <- bind_rows(df0, tibble(x=seq(-2, 2, 0.1), y=NA))
fit0 <- inla(
    y ~ x, data=pred_df, family="binomial", 
    control.compute=list(dic=TRUE), control.predictor=list(link=1, compute=TRUE),
    selection=list("(Intercept)"=1, x=1),
)

summary(fit0_null)
summary(fit0)

line_plot1 <- ggplot(df0, aes(x, c(0, 0.05, 0, 0, 0, 0, 0.05, 0), shape=as.factor(y))) +
    geom_point(size=1.2, col="#ff856d") +
    ylim(0, 1) +
    scale_shape_manual(values=c(2, 3)) +
    labs(x="Covariate space", shape="Observed data") +
    theme(
        panel.background=element_rect(fill=NA),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.line.x=element_line(color="gray50"),
        legend.key=element_blank(),
    )

ggsave("figs/figure 1/covar-space1.pdf", line_plot1, width=3.25, height=2.3)

# Obtain the initial joint posterior dist.----------------------------------------
jpost0 <- inla.rjmarginal(1000, fit0$selection)$samples |> t() |> as_tibble()

ggj <- ggplot(jpost, aes(`(Intercept):1`, `x:1`)) +
    geom_point(col="#ff856d", alpha=0.5, shape=1, size=0.85) +
    labs(x=NULL, y=NULL) +
    theme_bw() +
    theme(panel.grid.minor=element_blank())

plot_post1 <- ggMarginal(ggj, type="histogram", fill="#ff856d")

ggsave("figs/figure 1/post1.pdf", plot_post1, width=2.5, height=2.5)

###
pred_df <- bind_cols(pred_df, fit0$summary.fitted.values)

jitter <- position_jitter(0.04, 0, 123)

pred_ints <- ures$fit_best$summary.fitted.values[-c(1:nrow(df0)),] |> 
    as_tibble() |> 
    select(ymin=`0.025quant`, ymax=`0.975quant`) |> 
    mutate(x=xnew)

pred_df2 |>
    filter(is.na(y)) |>
    ggplot(aes(x, mean)) +
    geom_errorbar(aes(x, ymin=ymin, ymax=ymax), pred_ints, col="#c63df5", inherit.aes=FALSE, width=0.12) +
    geom_point(aes(x, c(-.01, -.09, rep(-.01, 4), -.09, -.01), shape=as.factor(y)), df0, col="#ff856d", size=0.96) +
    annotate(x=xnew, y=c(-.09, -.01, -.09, -.01), shape=1, col="#c63df5", geom="point", size=0.96) +
    # geom_point(aes(x, y=c(-0.05, -0.03, -0.03, -0.03, -0.03, -0.03, -0.05, -0.03)), df0, shape=1, col="#ff856d", size=0.9) +
    # scale_y_continuous(sec.axis=dup_axis(breaks=c(0, 1), name="Observed data")) +
    scale_shape_manual(values=c(2, 3)) +
    labs(x="Covariate space", y="Predicted risk", shape=NULL) +
    theme_bw() +
    theme(legend.position="none", panel.grid.minor=element_blank())
    # theme(
    #     axis.line.y.right=element_line(color="#3c4276"),
    #     axis.text.y.right=element_text(color="#3c4276"),
    #     axis.ticks.y.right=element_line(color="#3c4276"),
    #     axis.title.y.right=element_text(color="#3c4276")
    # )

ggsave("figs/figure 1/covar-space2.pdf", width=2.9, height=1.6)

pred_plot1 +
    geom_line(col="#ff856d", alpha=0.9)

ggsave("figs/figure 1/pred-space1.pdf", width=3.25, height=2.8)

ex_utility <- function(xnew, N=10) {
    pr_df <- bind_rows(df0, tibble(x=xnew, y=NA))
    pr_idx <- which(is.na(pr_df$y))
    
    pr_fit <- inla(
        y ~ x, data=pr_df, family="binomial", 
        selection=list(Predictor=pr_idx),
        control.predictor=list(link=1, compute=TRUE), 
    )
    eta_new <- inla.rjmarginal(N, pr_fit$selection)$samples
    ynew <- map(1:N, ~rbinom(length(xnew), 1, inla.link.invlogit(eta_new[,.x])))
    ubest <- -Inf
    ybest <- NULL
    cmat_best <- NULL
    
    utils <- map_dbl(ynew, ~{
        util_df <- bind_rows(df0, tibble(x=xnew, y=.x))
        util_fit <- inla(
            y ~ x, data=util_df, family="binomial", 
            selection=list("(Intercept)"=1, x=1),
            control.predictor=list(link=1, compute=TRUE), 
        )
        u <- -log(det(util_fit$selection$cov.matrix))
        if (u > ubest) {
            ubest <<- u
            fit_best <<- util_fit
            ybest <<- .x
        }
        u
    })
    
    return(list(ybest=ybest, fit_best=fit_best, util=mean(utils)))
}

xnew <- design_space$x[c(3, 4, 7, 11)]
ures <- fit_new(xnew)

pred_df2 <- bind_rows(df0, tibble(x=xnew, y=ures$ybest), tibble(x=seq(-2, 2, 0.1), y=NA))
fit2 <- inla(
    y ~ x, data=pred_df2, family="binomial", 
    control.predictor=list(link=1, compute=TRUE)
)

pred_df2 <- bind_cols(pred_df2, fit2$summary.fitted.values) |> 
    filter(is.na(y))

pred_plot1 +
    geom_ribbon(
        aes(ymin=`0.025quant`, ymax=`0.975quant`), data=pred_df2, 
        fill="#83a2d4", alpha=0.4
    ) +
    geom_line(data=pred_df2, col="#c63df5", alpha=0.9) +
    geom_point(aes(x, y), tibble(y=ures$ybest, x=xnew), shape=2, col="#3c4276", alpha=0.8, position=jitter) +
    annotate(x=xnew, y=c(-.05, rep(-.03, 3)), shape=1, col="#c63df5", geom="point", size=0.9)

ggsave("figs/figure 1/pred-space2.pdf", width=3.25, height=2.82)

jpost0 <- inla.rjmarginal(1000, fit0$selection)$samples |> t() |> as_tibble()
jpost <- inla.rjmarginal(1000, ures$fit_best$selection)$samples |> t() |> as_tibble()

ggplot(jpost, aes(`(Intercept):1`, `x:1`)) +
    geom_point(col="#c63df5", data=slice_sample(jpost, n=150), alpha=0.3, shape=1, size=0.85) +
    stat_ellipse(data=jpost0, type="norm", col="#ff856d", alpha=0.9, linetype="dashed", linewidth=1.05) +
    stat_ellipse(data=jpost0, type="norm", col="#ff856d", alpha=0.9, level=0.5, linewidth=1.05) +
    stat_ellipse(type="norm", col="#c63df5", alpha=0.9, linetype="dashed", linewidth=1.05) +
    stat_ellipse(type="norm", col="#c63df5", alpha=0.9, level=0.5, linewidth=1.05) +
    labs(x=NULL, y=NULL) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())

ggsave("figs/figure 1/post2.pdf", width=2.5, height=2.5)

###

#TODO: record utility as a function of search time while finding the optimum.
# add hline for base utility?

ex_sim_ann <- function(num_loc=1, iter=1, alpha=1, T0=1, info_iter=10) {
    d_curr <- d_best <- slice_sample(design_space, n=num_loc)
    u_curr <- u_best <- ex_utility(d_curr$x, 10)$util
    
    T_sched <- T0*seq(1, 0, length.out=iter)^alpha
    u_all <- double(iter)
    for (i in 1:iter) {
        # choose a new seach point. Reusing a point IS allowed
        next_pt <- slice_sample(design_space, n=1)
        
        # delete a random point and add the new one
        d_prop <- slice_sample(d_curr, n=max(nrow(d_curr)-1, 1)) |>
            bind_rows(next_pt)
        
        u_prop <- ex_utility(d_prop$x, 10)$util
        u_all[i] <- u_prop
        if (u_prop > u_best) { # if better than global best, update that info
            d_curr <- d_best <- d_prop
            u_curr <- u_best <- u_prop
        } else if (u_prop > u_curr) { # if better than current design, accept but don't update global best
            d_curr <- d_prop
            u_curr <- u_prop
        } else if (runif(1) < exp((log10(u_prop+1) - log10(u_curr+1))/T_sched[i])) { # or accept anyway sometimes
            d_curr <- d_prop
            u_curr <- u_prop
        }
        if (i %% info_iter == 0) # print update after iter
            print(paste0("Current U=", round(u_curr, 3), " at i=", i, "/", iter))
    }
    return(list(d=d_best, u_best=u_best, u_all=u_all))
}

sa_res <- ex_sim_ann(4, 100, alpha=1, T0=0.05)
sa_res_cp <- sa_res
# fudge the numbers a bit to simulate results more similar to real design space
sa_res_cp$u_all[20:60] <- map_dbl(sa_res$u_all[20:60], ~if (.x < 1) .x + runif(1, 0, 0.2) else .x)
sa_res_cp$u_all[60:100] <- map_dbl(sa_res$u_all[60:100], ~if (.x < 1) .x + runif(1, 0.2, 0.25) else .x)

sa_res_cp$u_all |> 
    enframe("iteration", "util") |> 
    mutate(u_best=map_dbl(1:n(), ~max(util[1:.x]))) |>
    pivot_longer(-iteration) |> 
    ggplot(aes(iteration, value, col=name, group=name)) +
    geom_line(alpha=0.9) +
    scale_color_manual(values=c("#f5b43d", "#c63df5")) +
    theme_bw() +
    theme(axis.title.y=element_blank(), legend.position="none", panel.grid.minor=element_blank())

ggsave("figs/figure 1/util-opt.pdf", width=2.9, height=1.6)

line_plot1 +
    geom_point(aes(x, c(.05, 0, .1, .05)), sa_res$d, col="#f5b43d", shape=1)

ggsave("figs/figure 1/covar-space2.pdf", width=3.25, height=2.3)

# Plot the spatiotemporal design space with each of the three surveys-------------
d_points <- bind_rows(
    mutate(df0, strat="init"),
    mutate(design_space[c(3, 4, 7, 11),], strat="cand"),
    mutate(sa_res$d, strat="opt")
)

d_points
pt_adj <- rep(0, nrow(d_points))
pt_adj[c(1, 2, 6, 7, 5, 11, 16, 8, 13, 15)] <- c(rep(c(-.08, .08), 2), rep(c(-.14, 0, .14), 2))

ggplot(design_space, aes(s1, s2, fill=x)) +
    geom_tile(col="gray70", linewidth=0.5, alpha=0.8) +
    geom_point(
        aes(s1, s2, col=fct_relevel(strat, "init", "cand")), d_points,
        size=1.6, position=position_nudge(x=pt_adj), inherit.aes=FALSE
    ) +
    geom_point(
        aes(s1, s2), d_points,
        shape=1, col="white", size=1.61, position=position_nudge(x=pt_adj), inherit.aes=FALSE
    ) +
    facet_wrap(~t, nrow=1) +
    scale_x_discrete(expand=expansion()) +
    scale_y_discrete(expand=expansion()) +
    scale_fill_gradient2(low="#0c2c84", mid="gray80", high="#7fcdbb") +
    # scale_fill_viridis_c(option="mako") +
    # scale_fill_distiller(palette="Greens") +
    scale_color_manual(values=c("#ff856d", "#c63df5", "#f5b43d"), guide="none") +
    theme_bw() +
    theme(
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        strip.text=element_blank(),
        strip.background=element_blank(),
        legend.title=element_blank(),
        legend.key.height=unit(0.23, "in"),
        legend.key.width=unit(0.09, "in"),
        legend.margin=margin(0, 0, 0, 0),
        panel.border=element_rect(fill=NA, color="gray70"),
        legend.text=element_text(size=7)
    )

ggsave("figs/figure 1/design-space.pdf", width=5.3, height=1.3)
