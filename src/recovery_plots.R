library(mcce)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(tidytext)
library(pmwg)
library(stringr)
library(paletti)
library(readr)

waterside_cols <- read_lines(file = "palette.txt")
viz_palette(waterside_cols, "waterside park")
fill_palette <- get_scale_fill(get_pal(waterside_cols))
col_palette <- get_scale_colour(get_pal(waterside_cols))

# Load all the samples
data_location <- here::here("data", "output", "5ModelRecovery")
recovery_files <- c(
  "NumericVDCE_CB_wvuhRfvyySjv_untagged.RData",
  "NumericVDCE_FPP_1878515.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_IEX_1878513.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_IST_1878369.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_MW_1878516.rcgbcm_5ModelRecovery.RData"
)

samples <- lapply(recovery_files, function(x) {
  print(paste("Extracting", x))
  get_samples(here::here(data_location, x))
})
names(samples) <- sapply(strsplit(recovery_files, "_"), "[[", 2)

original_samples <- get_samples(
  here::here(
    "data",
    "output",
    "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
  )
)

samples[["Original"]] <- original_samples

model_order <- c("CB", "FPP", "MW", "IST", "IEX", "Original")

recovery_data_files <- c(
  "NumericVDCE_CB_1878514.rcgbcm_5ModelRecovery_data.RDS",
  "NumericVDCE_FPP_1878515.rcgbcm_5ModelRecovery_data.RDS",
  "NumericVDCE_IEX_1878513.rcgbcm_5ModelRecovery_data.RDS",
  "NumericVDCE_IST_1878369.rcgbcm_5ModelRecovery_data.RDS",
  "NumericVDCE_MW_1878516.rcgbcm_5ModelRecovery_data.RDS"
)

recovery_data <- lapply(recovery_data_files, function(x) {
  print(paste("Reading data", x))
  readRDS(here::here(data_location, x))
})
names(recovery_data) <- sapply(strsplit(recovery_data_files, "_"), "[[", 2)

original_data <- readRDS(here::here("data", "output", "Task1_preprocessed.RDS"))
recovery_data[["Original"]] <- original_data
recovery_data <- bind_rows(recovery_data, .id = "source")

recovery_data %>%
  filter(rt < 5) %>%
  mutate(cell=paste0(price, rating)) %>%
  ggplot(aes(x = rt, colour = source)) +
  geom_density(size=1) +
  col_palette() +
  facet_wrap(~ cell)

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
  filter(source == "Original", Parameter == "IEX") %>%
  arrange(desc(rel_val)) %>%
  pull(subjectid)

Par_order <- model_medians %>%
  filter(source == "Original") %>%
  group_by(Parameter) %>%
  summarise(median_val = median(rel_val)) %>%
  arrange(median_val) %>%
  pull(Parameter)

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
  mutate(subjectid = factor(as.numeric(subjectid), levels = 1:27, labels = c(letters, "Group"))) %>%
  mutate(highlight = ifelse(subjectid == "Group", "Highlight", "Normal")) %>%
  mutate(Parameter = factor(Parameter, Par_order)) %>%
  ggplot(aes(x = subjectid, y = rel_val, fill = Parameter, alpha = highlight)) +
    geom_col() +
    xlab("Subject Identifier") +
    ylab("Relative Evidence") +
    fill_palette() +
    scale_alpha_manual(values=c(1,0.8)) +
    facet_grid(rows = vars(source)) +
    scale_y_continuous(labels = NULL, breaks = NULL)

model_medians %>%
  mutate(subjectid = factor(subjectid, subject_order)) %>%
  mutate(subjectid = fct_relevel(subjectid, "Group", after=Inf)) %>%
  mutate(subjectid = group_indices(subjectid))

  mutate(Parameter = factor(Parameter, Par_order)) %>%
  ggplot(aes(x = Parameter, y = rel_val, fill = Parameter)) +
  geom_col() +
  fill_palette() +
  labs(y = "Relative Evidence") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(rows = vars(source), cols = vars(subjectid))


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
    geom_point(size = 0.5, colour = waterside_cols[3]) +
    geom_point(data = recovered %>% filter(subjectid == "Group"), colour = waterside_cols[6], size = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(vars(Parameter), scales = "free") +
    labs(
      x = "Generating Value",
      y = "Recovered Value",
      title = paste("Generated from", model),
      caption = scatter_caption
    ) +
    scatter_theme +
    col_palette()
    print(p)
}
