library(tidyverse)
library(mcce)
library(patchwork)
library(pmwg)
library(watercolours)

# original bigrun
filenames <- list.files(
  here::here("data", "output", "bigrun"),
  pattern = ".RData",
  full.names = TRUE
)
samplers <- sapply(filenames, get_samples, simplify = FALSE)
names(samplers) <- substring(names(samplers), 114, 117)

menu_args <- list(
  choices = c(samplers[[1]]$par_names, "Quit"),
  title = "What parameter do you want a traceplot of?"
)

# Traceplots of individual theta_mu parameters across the 3 runs
tmus <- sapply(names(samplers), function(x) {
                 y <- samplers[[x]]
                 y %>%
                   extract_tmu %>%
                   pivot_longer(cols = -c(sampleid, stageid),
                                names_to = "par") %>%
                   rename_with(.cols = c(stageid, value),
                               .fn = function(z) {
                                 paste0(z, ".", x)
                                })},
               simplify = FALSE) %>%
  purrr::reduce(left_join, by = c("sampleid", "par")) %>%
  pivot_longer(contains("run"),
               names_to = c(".value", "run"),
               names_pattern = "(.*).run(.*)") %>%
  mutate(run = factor(run),
         stageid = factor(stageid, levels = c("sample", "burn", "init", "adapt")))

while ((menu_choice <- do.call(menu, menu_args)) != 22) {
  par_name <- samplers[[1]]$par_names[menu_choice]

  par_samples <- tmus %>%
    filter(par == par_name)

  transitions <- par_samples %>%
    filter(stageid != lag(stageid)) %>%
    select(sampleid)

  g <- par_samples %>%
    ggplot(aes(x = sampleid, y = value, colour = run, linetype = stageid)) +
    geom_line() +
    scale_colour_watercolour() +
    labs(title = paste("Parameter:", par_name)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank()) +
    geom_vline(data = transitions, mapping = aes(xintercept = sampleid),
               lty = "dotted")
    print(g)
}

cor_matrices <- sapply(samplers, FUN = function(x) {
  extract_cov(x, filter = "sample") %>%
    apply(c(1,2), mean) %>%
    cov2cor
}, simplify = FALSE)

cor_df <- sapply(names(cor_matrices), FUN = function(x) {
  cor_matrices[[x]] %>%
    matrix_to_df %>%
    mutate(run = x)
}, simplify = FALSE) %>%
  bind_rows()


cor_df %>%
  matdf_plot() +
  facet_wrap(~ run) +
  scale_fill_watercolour(type = "diverging")

## Differences
cor_df %>%
  group_by(X, Y) %>%
  mutate(mean_corr = mean(value)) %>%
  ungroup %>%
  mutate(value = abs(value - mean_corr)) %>%
  matdf_plot() +
  facet_wrap(~ run) +
  scale_fill_watercolour(type = "continuous") +
  ggtitle("Differences between each runs correlation to mean of all runs correlation")

tmus %>%
  filter(str_detect(par, "^alpha"), stageid == "sample") %>%
  group_by(par, run) %>%
  summarise(median_val = median(value)) %>%
  pivot_wider(names_from=par, values_from=median_val)

pairwise <- list()
pairwise[["run1 - run2"]] <- cor_matrices[[1]] - cor_matrices[[2]]
pairwise[["run2 - run3"]] <- cor_matrices[[2]] - cor_matrices[[3]]
pairwise[["run1 - run3"]] <- cor_matrices[[1]] - cor_matrices[[3]]

pairwise_df <- sapply(names(pairwise), FUN = function(x) {
  pairwise[[x]] %>%
    matrix_to_df %>%
    mutate(comparison = x)
}, simplify = FALSE) %>%
  bind_rows()

pairwise_df %>%
  matdf_plot() +
  facet_wrap(~ comparison) +
  scale_fill_watercolour(type = "diverging")
