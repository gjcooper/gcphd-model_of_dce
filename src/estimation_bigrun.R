library(tidyverse)
library(mcce)
library(patchwork)
library(pmwg)
library(watercolours)

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
while ((menu_choice <- do.call(menu, menu_args)) != 22) {
  par_name <- samplers[[1]]$par_names[menu_choice]

  par_samples <- sapply(names(samplers), FUN = function(x) {
    samplers[[x]] %>%
    as_mcmc %>%
    data.frame %>%
    tibble %>%
    mutate(sample_id = row_number()) %>%
    mutate(stage = samplers[[x]]$samples$stage) %>%
    pivot_longer(cols = -c(sample_id, stage), names_to = "parameter") %>%
    filter(parameter == par_name) %>%
    rename_with(.cols = c(stage, value), .fn = function(y) {paste0(y, ".", x)})
  }, simplify = FALSE)

  par_samples <- purrr::reduce(
    par_samples,
    left_join,
    by = c("sample_id", "parameter")
  ) %>%
  pivot_longer(
    contains("run"),
    names_to = c(".value", "run"),
    names_pattern = "(.*).run(.*)"
  ) %>%
  mutate(run = factor(run)) %>%
  mutate(stage = factor(stage, levels = c("sample", "burn", "init", "adapt")))

  transitions <- par_samples[par_samples$stage != lag(par_samples$stage), "sample_id"]
  transitions <- transitions[complete.cases(transitions),]

  g <- par_samples %>%
    ggplot(aes(x = sample_id, y = value, colour = run, linetype = factor(stage))) +
    geom_line() +
    scale_colour_watercolour() +
    labs(title = paste("Parameter:", par_name)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank()) +
    geom_vline(data = transitions, mapping = aes(xintercept = sample_id), lty = "dotted")
    print(g)
}


cov_matrices <- sapply(names(samplers), FUN = function(x) {
  samplers[[x]] %>%
  as_mcmc(selection = "theta_sig", filter = "sample") %>%
  data.frame %>%
  tibble %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(everything(), names_to = c("X", "Y"), names_sep = "\\.") %>%
  mutate(across(c(X, Y), factor)) %>%
  xtabs(value ~ X + Y, .)
}, simplify = FALSE)

cor_matrices <- sapply(names(cov_matrices), FUN = function(x) {
  cov_matrices[[x]] %>%
  cov2cor()
}, simplify = FALSE)

cor_df <- sapply(names(cor_matrices), FUN = function(x) {
  cor_matrices[[x]] %>%
    data.frame %>%
    tibble %>%
    rename(!!x := Freq)
}, simplify = FALSE) %>%
  purrr::reduce(left_join, by = c("X", "Y")) %>%
  pivot_longer(c(run1, run2, run3), names_to = "run", values_to = "correlation")

cor_df %>%
  ggplot(aes(X, Y, fill = correlation, label=round(correlation,2))) +
  geom_tile() +
  geom_text(size=2) +
  facet_wrap(~ run) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_equal() +
  scale_fill_watercolour(discrete=FALSE)

## Differences
cor_df %>%
  group_by(X, Y) %>%
  mutate(mean_corr = mean(correlation)) %>%
  ungroup %>%
  mutate(corr_diff = abs(correlation - mean_corr)) %>%
  ggplot(aes(X, Y, fill = corr_diff)) +
  geom_tile() +
  facet_wrap(~ run) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_equal() +
  scale_fill_watercolour(type = "continuous") +
  ggtitle("Differences between each runs correlation to mean of all runs correlation")

arch_samples <- sapply(names(samplers), FUN = function(x) {
  samplers[[x]] %>%
    as_mcmc(filter = "sample") %>%
    data.frame %>%
    tibble %>%
    select(starts_with("alpha")) %>%
    summarise(across(everything(), median))
}, simplify=FALSE) %>% print


pairwise <- list()
pairwise[["run1 - run2"]] <- cor_matrices[[1]] - cor_matrices[[2]]
pairwise[["run2 - run3"]] <- cor_matrices[[2]] - cor_matrices[[3]]
pairwise[["run1 - run3"]] <- cor_matrices[[1]] - cor_matrices[[3]]

pairwise_df <- sapply(names(pairwise), FUN = function(x) {
  pairwise[[x]] %>%
    data.frame %>%
    tibble %>%
    rename(!!x := Freq)
}, simplify = FALSE) %>%
  purrr::reduce(left_join, by = c("X", "Y")) %>%
  pivot_longer(c(`run1 - run2`, `run2 - run3`, `run1 - run3`), names_to = "comparison", values_to = "difference")

pairwise_df %>%
  ggplot(aes(X, Y, fill = difference, label=round(difference,2))) +
  geom_tile() +
  geom_text(size=2) +
  facet_wrap(~ comparison) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_equal() +
  scale_fill_watercolour(type = "diverging")
