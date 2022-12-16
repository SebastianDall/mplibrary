library(ggplot2)

# test <- stat.test
stages_all <- c("Pre", "treatment_5", "treatment_10", "treatment_15", "treatment_21", "Post", "followup_1m", "followup_3m", "followup_6m", "followup_12m", "Donor")



#' Similiarity to donor plot
#' #'
#' Creates a similarity to donor plot.
#'
#' @param metadata_beta_summarised A summarised beta diversity df
#' @param method Which method was used (either sorensen or bray)
#'
#' @export
#'
#'
gg_beta <- function(metadata_beta_summarised, method) {
    if (method == "bray") {
        y_lab <- "Similarity to Donors (Bray-Curtis)"
    } else if (method == "sorensen") {
        y_lab <- "Similarity to Donors (Sørensen Coefficient)"
    } else {
        stop("method is not valid.")
    }

    gg <- metadata_beta_summarised %>%
        ggplot(aes(x = x_axis, y = 1 - median_dissimilarity, color = x_axis)) +
        geom_boxplot(outlier.shape = NA) +
        geom_line(aes(group = id), alpha = 0.6, color = "grey70", position = position_dodge(0.2)) +
        geom_point(aes(group = id), position = position_dodge(0.2), size = 2) +
        labs(y = y_lab, x = "") +
        facet_grid(. ~ paste0(group, " - ", actual_donor)) +
        plot_theme

    return(gg)
}
