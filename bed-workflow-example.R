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

ggplot(design_space, aes(s1, s2, fill=x)) +
    geom_tile(col="gray70", linewidth=0.6) +
    geom_point(
        aes(s1, s2), df0,
        shape=1, col="#ff856d", size=1.1, position=position_nudge(x=c(-0.1, 0.1, 0, 0, 0, 0.1, -0.1, 0)), inherit.aes=FALSE
    ) +
    facet_wrap(~t, nrow=1) +
    scale_x_discrete(expand=expansion()) +
    scale_y_discrete(expand=expansion()) +
    # scale_fill_viridis_c(option="mako") +
    scale_fill_distiller(palette="Greens") +
    theme_bw() +
    theme(
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        strip.text=element_blank(),
        strip.background=element_blank(),
        legend.title=element_blank(),
        # legend.key.height=unit(1, "in"),
        legend.key.width=unit(0.1, "in"),
        legend.margin=margin(0, 0, 0, 0),
        panel.border=element_rect(fill=NA, color="gray70")
        # legend.text=element_text(size=4)
    )

ggsave("figs/figure 1/design-space.pdf", width=5.3, height=1.5)

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
    

pred_df <- bind_cols(pred_df, fit0$summary.fitted.values)

jitter <- position_jitter(0.04, 0, 123)

pred_plot1 <- pred_df |>
    filter(is.na(y)) |>
    ggplot(aes(x, mean)) +
    geom_ribbon(aes(ymin=`0.025quant`, ymax=`0.975quant`), fill="#9ac9e7", alpha=0.4) +
    geom_point(aes(y=y), df0, shape=2, col="#3c4276", alpha=0.8, position=jitter) +
    geom_point(aes(x, y=c(-0.05, -0.03, -0.03, -0.03, -0.03, -0.03, -0.05, -0.03)), df0, shape=1, col="#ff856d", size=0.9) +
    scale_y_continuous(sec.axis=dup_axis(breaks=c(0, 1), name="Observed data")) +
    labs(x="Covariate space", y="Predicted risk") +
    theme_bw() +
    theme(
        axis.line.y.right=element_line(color="#3c4276"),
        axis.text.y.right=element_text(color="#3c4276"),
        axis.ticks.y.right=element_line(color="#3c4276"),
        axis.title.y.right=element_text(color="#3c4276")
    )

pred_plot1 +
    geom_line(col="#ff856d", alpha=0.9)

ggsave("figs/figure 1/pred-space1.pdf", width=3.25, height=2.8)

fit_new <- function(xnew) {
    pr_df <- bind_rows(df0, tibble(x=xnew, y=NA))
    pr_idx <- which(is.na(pr_df$y))
    
    pr_fit <- inla(
        y ~ x, data=pr_df, family="binomial", 
        selection=list(Predictor=pr_idx),
        control.predictor=list(link=1, compute=TRUE), 
    )
    eta_new <- inla.rjmarginal(10, pr_fit$selection)$samples
    ynew <- map(1:10, ~rbinom(length(xnew), 1, inla.link.invlogit(eta_new[,.x])))
    ubest <- -Inf
    ybest <- NULL
    fit_best <- NULL
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
            cmat_best <<- util_fit$selection$cov.matrix
            fit_best <<- util_fit
            ybest <<- .x
        }
        u
    })
    
    return(list(ybest=ybest, fit_best=fit_best, cmat_best=cmat_best, util=mean(utils)))
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
    geom_point(col="gray50", alpha=0.5, shape=1) +
    stat_ellipse(data=jpost0, type="norm", col="#ff856d", alpha=0.8, linewidth=1.05) +
    stat_ellipse(type="norm", col="#c63df5", alpha=0.8, linewidth=1.05) +
    labs(x=NULL, y=NULL) +
    theme_bw()
