library(mcce)
library(dplyr)
library(ggplot2)
#library(tidyr)
#library(forcats)
#library(pmwg)
library(stringr)
library(patchwork)

# Load all the samples
accept_file <- here::here("data", "output", "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData")
reject_file <- here::here("data", "output", "NumericVDCE_2117450.rcgbcm_RejectEstimation5Model.RData")

accept_samples <- get_samples(accept_file)
reject_samples <- get_samples(reject_file)

samples <- list(accept = accept_samples, reject = reject_samples)

model_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha")) %>%
    get_medians() %>%
    group_by(subjectid) %>%
    mutate(rel_val = value / sum(value))
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

model_medians <- bind_rows(model_medians, .id = "source") %>%
  mutate(source = factor(source, levels = c("accept", "reject"))) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

model_plot <- function(medians) {
  medians %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    xlab("Subject Identifier") +
    ylab("Relative Evidence") +
    scale_fill_brewer(palette = "Dark2", name = "Model") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

accept_plot <- model_medians %>% filter(source == "accept")  %>% model_plot
reject_plot <- model_medians %>% filter(source == "reject")  %>% model_plot
accept_plot / reject_plot

par_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha", negate = TRUE)) %>%
    get_medians(alpha = FALSE)
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

par_medians <- bind_rows(par_medians, .id = "source") %>% mutate(value = log(value))

par_medians %>%
  ggplot(aes(x = Parameter, y = exp(value), color = source)) +
  geom_boxplot()

accept_pars <- par_medians %>% filter(source == "accept")
reject_pars <- par_medians %>% filter(source == "reject")

combined <- recovery %>%
  left_join(original, by = c("Parameter", "subjectid")) %>%
  select(-source.y) %>%
  rename(
    recovered_value = value.x,
    estimated_value = value.y,
    recovery_model = source.x
  ) %>%
  mutate(source = factor(recovery_model, levels = c("accept", "reject"))) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

scatter_theme <- theme(
  plot.title = element_text(hjust = 0.5),
  plot.caption = element_text(hjust = 0),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
)
scatter_caption <- str_wrap(paste(
  "The data generating values are the medians of the posterior for the full",
  "model estimated with the dirichlet process. The generation process uses the",
  "random sample function from only one of the 5 models.",
  "The recovered values are again the medians of the posterior for the full",
  "model (again with the dirichlet process)."
), 100)
scatter_caption <- paste(
  scatter_caption,
  "\n",
  str_wrap(paste(
    "Shown in red is the theta_mu value from the original data and the",
    "recovered data, but is not actually used to generate anything"
    ), 100)
)

for (model in c("accept", "reject")) {
  recovered <- combined %>% filter(recovery_model == model)
  p <- ggplot(recovered, aes(x = estimated_value, y = recovered_value)) +
    geom_point(size = 0.5) +
    geom_point(data = recovered %>% filter(subjectid == "Group"), colour = "red", size = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(vars(Parameter), scales = "free") +
    labs(
      x = "Generating Value",
      y = "Recovered Value",
      title = paste("Generated from", model),
      caption = scatter_caption
    ) +
    scatter_theme
    print(p)
}

