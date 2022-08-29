library(mcce)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(tidytext)
library(pmwg)
library(stringr)
library(watercolours)
library(readr)
library(patchwork)

task <- "Preferential" ## Other option, something like Preferential

if (task == "Veridical") {
  # Load all the samples
  data_location <- here::here("data", "output", "VeridicalRecovery")
  #data_location <- here::here("data", "output", "VeridicalRecovery_Jun2021")
  odata <- "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
} else {
  # Load all the samples
  data_location <- here::here("data", "output", "PrefRecovery")
  #data_location <- here::here("data", "output", "PrefRecovery_Oct2021")
  odata <- "PrefDCE_2506730.rcgbcm_Estimation5Model.RData"
}
recovery_files <- list.files(data_location, pattern = ".RData")
recovery_data_files <- list.files(data_location, pattern = ".RDS")

samples <- lapply(recovery_files, function(x) {
  print(paste("Extracting", x))
  get_samples(here::here(data_location, x))
})
names(samples) <- sapply(strsplit(recovery_files, "_"), "[[", 2)

original_samples <- get_samples(
  here::here(
    "data",
    "output",
    odata
  )
)

samples[["Original"]] <- original_samples

model_order <- c("CB", "FPP", "MW", "IST", "IEX", "Original")

recovery_data <- lapply(recovery_data_files, function(x) {
  print(paste("Reading data", x))
  readRDS(here::here(data_location, x))
})
names(recovery_data) <- sapply(strsplit(recovery_data_files, "_"), "[[", 2)

original_data <- original_samples$data
recovery_data[["Original"]] <- original_data
recovery_data <- bind_rows(recovery_data, .id = "source")

recovery_data %>%
  filter(rt < 10) %>%
  mutate(accept = factor(accept, labels = c("Reject", "Accept"))) %>%
  ggplot(aes(x = source, y = rt, colour = generator)) +
  geom_boxplot() +
  facet_wrap(~ accept) +
  scale_colour_manual(name = "Generator Model", values = unname(frankwebb_cols))


recovery_data %>%
  filter(rt < 5) %>%
  mutate(simulated = ifelse(source == "Original", "No", "Yes")) %>%
  mutate(cell=paste0(price, rating)) %>%
  ggplot(aes(x = rt, colour = source, linetype=simulated)) +
  geom_density(size=0.7) +
  frank_colmap2 +
  facet_wrap(~ cell)

ggsave(
  filename = here::here(
    "results",
    task,
    paste0("SimulatedRT_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

model_medians <- sapply(samples, function(x) {
    extract_parameters(x, str_subset(x$par_names, "alpha")) %>%
      get_medians() %>%
      group_by(subjectid) %>%
      mutate(rel_val = value / sum(value))
    },
    USE.NAMES = TRUE,
    simplify = FALSE
  ) %>%
  bind_rows(.id = "source") %>%
  mutate(source = factor(source, levels = model_order)) %>%
  mutate(Parameter = as.character(Parameter)) %>%
  mutate(subjectid = case_when(
    subjectid == "theta_mu" ~ "Group",
    TRUE ~ str_pad(subjectid, 2, pad = "0")
  ))

subject_order <- model_medians %>%
  filter(subjectid != "Group") %>%
  group_by(subjectid, source) %>%
  summarise(maxRel = max(rel_val)) %>%
  group_by(subjectid) %>%
  summarise(meanMax = mean(maxRel)) %>%
  arrange(desc(meanMax)) %>%
  pull(subjectid)

Par_order <- model_medians %>%
  filter(source == "Original") %>%
  group_by(Parameter) %>%
  summarise(median_val = median(rel_val)) %>%
  arrange(median_val) %>%
  pull(Parameter)

subject_recovery <- model_medians %>%
  filter(source != "Original") %>%
  filter(subjectid != "Group") %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  mutate(source = factor(source, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    labs(y = "Relative Evidence") +
    frank_colmap +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank()) +
    facet_grid(rows = vars(source)) +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    ggtitle("Indidividual Subjects")

group_recovery <- model_medians %>%
  filter(source != "Original") %>%
  filter(subjectid == "Group") %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  mutate(source = factor(source, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    frank_colmap +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    facet_grid(rows = vars(source)) +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    ggtitle("Group Level")

subject_recovery + group_recovery +
  plot_layout(widths = c(6, 1), guides = "collect") +
  plot_annotation(
    title = "Relative Evidence for each of the 5 architectures",
    subtitle = "Each row corresponds to the generating architecture",
    caption = "The area of each stacked bar is the relative proportion of each of the five possible modelled architectures as their dirichlet process"
    )

ggsave(
  filename = here::here(
    "results",
    task,
    paste0("ArchRecovery_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

par_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha", negate = TRUE)) %>%
    get_medians(alpha = FALSE)
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

par_medians <- bind_rows(par_medians, .id = "source") %>% mutate(value = log(value))

recovery <- par_medians %>% filter(source != "Original")
original <- par_medians %>% filter(source == "Original")

combined <- recovery %>%
  left_join(original, by = c("Parameter", "subjectid")) %>%
  select(-source.y) %>%
  rename(
    recovered_value = value.x,
    estimated_value = value.y,
    recovery_model = source.x
  ) %>%
  mutate(source = factor(recovery_model, levels = model_order[-6])) %>%
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

for (model in model_order[-6]) {
  recovered <- combined %>% filter(recovery_model == model)
  p <- ggplot(recovered, aes(x = estimated_value, y = recovered_value)) +
    geom_point(size = 0.5, colour = frankwebb_cols[3]) +
    geom_point(data = recovered %>% filter(subjectid == "Group"), colour = frankwebb_cols[6], size = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(vars(Parameter), scales = "free") +
    labs(
      x = "Generating Value",
      y = "Recovered Value",
      title = paste("Generated from", model),
      caption = scatter_caption
    ) +
    scatter_theme +
    frank_colmap2
    print(p)

  ggsave(
    filename = here::here(
      "results",
      task,
      paste0("ParRecovery_", model, "_", Sys.Date(), ".png")
    ),
    dpi = 200,
    width = 14.1,
    height = 7.53,
    units = "in"
  )
}
