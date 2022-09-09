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
