library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(paletti)
library(readr)

frankwebb_cols <- read_lines(file = "palette.txt")
viz_palette(frankwebb_cols, "Frank Webb palette")
names(frankwebb_cols) <- names_ll()
frank_colmap <- scale_fill_manual(name = "Architecture", values = frankwebb_cols)

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
  Par_order <- c("IST", "CB", "FPP", "MW", "IEX")
  medians %>%
    filter(subjectid != "Group") %>%
    mutate(Parameter = factor(Parameter, Par_order)) %>%
    ggplot(aes(x = subjectid, y = rel_val, fill = Parameter)) +
    geom_col() +
    frank_colmap +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Relative Evidence") +
    scale_y_continuous(labels = NULL, breaks = NULL)
}

design <- "
  11111
  22333
"
accept_plot <- model_medians %>% filter(source == "accept")  %>% model_plot
reject_plot <- model_medians %>% filter(source == "reject")  %>% model_plot
accept_plot + reject_plot + guide_area() + plot_layout(design = design, guides = "collect")

ggsave(
  filename = here::here(
    "results",
    "Veridical",
    paste0("EstimatedArch_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)


par_medians <- sapply(samples, function(x) {
  extract_parameters(x, str_subset(x$par_names, "alpha", negate = TRUE)) %>%
    filter(subjectid != "theta_mu") %>%
    get_medians(alpha = FALSE)
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

par_medians <- bind_rows(par_medians, .id = "source") %>% mutate(value = log(value))

par_colours <- c("grey", frankwebb_cols[3], frankwebb_cols[2], "grey", frankwebb_cols[3], frankwebb_cols[3], frankwebb_cols[2], frankwebb_cols[2], frankwebb_cols[3], frankwebb_cols[3], frankwebb_cols[2], frankwebb_cols[2], frankwebb_cols[3], frankwebb_cols[3], frankwebb_cols[2], frankwebb_cols[2])

gencols <- scale_colour_manual(name = "Generating Architecture", values = frankwebb_cols[c(2, 1)])

par_medians %>%
  mutate(source = recode(source, accept = "IEX", reject = "IST")) %>%
  mutate(colour = rep(par_colours, n_distinct(.$subjectid))) %>%
  ggplot(aes(x = Parameter, y = exp(value), fill = colour, color = source)) +
  gencols +
  scale_fill_identity() +
  geom_boxplot()

ggsave(
  filename = here::here(
    "results",
    "Veridical",
    paste0("EstimatedPars_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)
