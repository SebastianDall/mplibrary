library(ggplot2)

# p.label = "p.format"
# p.list = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p < 0.0001", "p < 0.001", " p < 0.01", "p < 0.05", "ns"))

# test <- stat.test
stages_all <- c("Pre", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "Post", "followup_1m", "followup_3m", "followup_6m", "followup_12m", "Donor")


#' Title
#' #'
#' Description
#'
#' @param df A dataframe containing the columns x, alpha-diversity metric
#' @param x A vector with the variables from the x_axis column which will be plotted on the x-axis. The vector must be given in the order to appear on the x-axis.
#' @param alpha_metric A variable stating the name of the column, which contains the alpha-diversity metric.
#' @export
#'
#' @example
#'
gg_alpha <- function(df, x, alpha_metric = "richness") { # , stat.test = "wilcox.test"
    stopifnot(!is.null(x), "the stages to be plotted on the x-axis must be stated!")

    df_rename <- df %>%
        rename(alpha_metric = alpha_metric)

    gg_richness_data <- df_rename %>%
        filter(!is.na(remission) | x_axis == "Donor", x_axis %in% x) %>% #
        mutate(x_axis = factor(x_axis, levels = x)) %>%
        arrange(id)


    gg_alpha_rich <- gg_richness_data %>%
        ggplot(aes(x = x_axis, y = alpha_metric, color = x_axis)) +
        geom_point(aes(group = id), position = position_dodge(0.2)) +
        geom_boxplot(outlier.shape = NA) +
        geom_line(data = filter(gg_richness_data, x_axis != "Donor"), aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) +
        geom_point(aes(group = id), position = position_dodge(0.2), size = 2)

    return(gg_alpha_rich)
}

# stat_compare_means(label = p.label, method = test, comparisons = list(c("Pre", "Post")), paired = T, label.y = label.y, symnum.arg = p.list) +
# stat_compare_means(data = filter(gg_richness_data, group != "placebo"), label = p.label, method = test, comparisons = list(c("Pre", "Donor"), c("Post", "Donor")), symnum.arg = p.list)
# coord_cartesian(ylim = ylim)
# facet_grid(. ~ group, scales = "free_x", space = "free")

#   } else {


#     gg_richness_data <- df_rename %>%
#       mutate(x_axis = factor(x_axis, levels = stages_all, labels = stages_labels)) %>%
#       mutate(lty = ifelse(x_axis == "Donor", NA, "solid"),
#              shape = ifelse(x_axis == "Donor", NA, "pt_sample")) %>%
#       arrange(id)

#     gg_alpha_rich <- gg_richness_data %>%
#       ggplot(aes(x = x_axis, y = alpha_metric, color = id)) +
#       geom_line(aes(group = id, linetype = lty), alpha = 0.6, color = "grey70", show.legend = F) +
#       geom_point(aes(group = id, shape = shape), size = 2, show.legend = F) +
#       geom_point(data = filter(gg_richness_data, x_axis == "Donor"), aes(group = id), size = 2, alpha = 0.6, position = position_jitter()) +
#       coord_cartesian(ylim = ylim) +
#       plot_theme +
#       theme(legend.position = "right",
#             axis.text.x = element_text(angle = 90, vjust = 0.35, hjust = 1))


#     }
