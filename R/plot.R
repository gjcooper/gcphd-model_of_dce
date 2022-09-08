#' Plot subject and group level summaries of estimates
#'
#' This function takes subject level and group level summary statistics
#' from an mcce sampling run and creates a plot showing the results.
#'
#' @param alpha_summ A data.frame containing subject level summary stats of
#'   each parameter
#' @param tmu_summ A data.frame containing group level summary stats
#' @param transform An optional function to transform the summary stats.
#'
#' @return A ggplot object
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @export
plot_summary <- function(alpha_summ, tmu_summ, transform=identity) {
  median_alpha %>%
    mutate(across(-SubjectID, transform)) %>%
    pivot_longer(!SubjectID) %>%
    ggplot(mapping = aes(x = name, y = value)) +
    geom_point(colour = "#D0781C") +
    geom_point(
      data = pivot_longer(median_theta_mu %>% transform, everything()),
      mapping = aes(x = name, y = value),
      size = 3,
      colour = "black"
    ) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0))
}
