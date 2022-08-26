plot_medians <- function(median_alpha, median_theta_mu, transform=identity) {
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
