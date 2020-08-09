# Set of functions that have been used in prior work. Mostly semi-customized R plotting scripts. May be useful for 
# future work. 

#' Useful for Bayesian binomial models in which there are random effects by some grouping variable. Assumes an 
#' interest in plotting the raw data by group as points, with sizes varying as a function of the number of observations
#' per group. Creates a 95% CI around linear fit line using the stat_lineribbon function. Other features of the 
#' function may be overly customized for initial use case in summarizing DeYoung thesis data. But can be made more 
#' general if need be. 
create_single_percent_plot <- function(predict_df, data_path, fit, yvar, ylim, y_interval, y_hit, xvar, xlim, 
                                       x_interval, grp_var, grp_n_min, draws, color_pal, plot_title, x_label = NULL, 
                                       out_path = NULL) {
    require(tidyverse)
    require(glue)
    require(modelr)
    require(tidybayes)
    
    load(data_path)
    y_breaks <- seq(ylim[1], ylim[2], by = y_interval)
    y_labels <- paste0(y_breaks, "%")
    
    if(!is.null(x_label)){
        x_breaks <- seq(xlim[1], xlim[2], by = x_interval)
        x_labels <- seq(floor(xlim[1]), floor(xlim[2]), by = x_interval)
    }
    
    
    
    # Expand and create a prediction data frame to pull 
    model_pred <- predict_df %>% 
        add_fitted_draws(fit, n=draws, allow_new_levels = TRUE) %>% 
        ungroup() %>% 
        mutate(
            perc = .value*100, 
        ) %>% 
        select(!!sym(xvar), perc)
    
    # because dplyr is dumb when it comes to strings...
    model_pred[[xvar]] <- model_pred[[xvar]] + xlim[1]
    
    # Create a raw data set of outputs set minimum group N if NULL defaults to 0
    if(is.null(grp_n_min)){
        grp_n_min <- 0
    }
    
    out_by_grp <- dat_grads_model %>% 
        select(!!sym(xvar), !!sym(grp_var), !!sym(yvar)) %>% 
        group_by(!!sym(grp_var), !!sym(xvar)) %>% 
        summarize(
            perc = sum(!!sym(yvar) == y_hit, na.rm=TRUE)/sum(!is.na(!!sym(yvar)))*100, 
            count = sum(!is.na(!!sym(yvar)))
        ) %>% 
        filter(count > grp_n_min)
    
    # dplyr remains dumb when it comes to strings... 
    out_by_grp[[xvar]] <- out_by_grp[[xvar]] + xlim[1]
    
    g <- ggplot(model_pred, aes(x = !!sym(xvar), y = perc)) +
        stat_lineribbon(.width = .95, fill=RColorBrewer::brewer.pal(9, color_pal)[2], alpha = .85,
                        point_interval = mean_qi, color = RColorBrewer::brewer.pal(9, color_pal)[9]) +
        geom_jitter(data=out_by_grp, aes(x = !!sym(xvar), y = perc, size = count), alpha = .1875, width = .075, 
                    color = RColorBrewer::brewer.pal(9, color_pal)[9], show.legend = FALSE) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels) +
        coord_cartesian(ylim = ylim) +
        ggtitle(plot_title) + 
        theme_bw() +
        theme(axis.title.y = element_blank(), 
              plot.title = element_text(face = "italic"))
    
    if(!is.null(x_label)) {
        g <- g + scale_x_continuous(breaks = x_breaks, labels = x_labels) +
            labs(x = x_label)
    } else { 
        g <- g + theme(axis.ticks.x = element_blank(), 
                       axis.title.x = element_blank(), 
                       axis.text.x = element_blank())
    }
    
    if(!is.null(out_path)) {
        ggsave(out_path, plot = g, device = "png", dpi = 900, width = 9, height = 6, units = "in")
    }
    
    return(g)
}
