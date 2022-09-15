#' Plot subject and group level summaries of estimates
#'
#' This function takes subject level and group level summary statistics
#' from an mcce sampling run and creates a plot showing the results.
#'
#' @param summ A data.frame containing subject and group level summary stats
#'   for each parameter. The data.frame should have a subjectid column with
#'   theta_mu representing the group level estimates.
#' @param transform An optional function to transform the summary stats.
#'
#' @return A ggplot object
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @export
plot_summary <- function(summ, transform=identity) {
  tsumm <- summ %>%
    mutate(across(-subjectid, transform)) %>%
    pivot_longer(!subjectid)

  tsumm %>%
    filter(subjectid != "theta_mu") %>%
    ggplot(mapping = aes(x = name, y = value)) +
    geom_point(colour = "#D0781C") +
    geom_point(
      data = tsumm %>% filter(subjectid == "theta_mu"),
      mapping = aes(x = name, y = value),
      size = 3,
      colour = "black"
    ) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0))
}

#' Plot matrix_to_df object using ggplot
#'
#' This function takes a matrix in data.frame form, usually output from
#' matrix_to_df and displays it using ggplot.
#' The resulting figure can be adjusted further.
#'
#' @param x A data.frame containing an X and Y column for parameter labels
#'   and a value column containing the data to plot.
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @export
matdf_plot <- function(x) {
  ggplot(x, aes(X, Y, fill = value, label = round(value, 2))) +
    geom_tile() +
    coord_equal() +
    geom_text(size = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank())
}
