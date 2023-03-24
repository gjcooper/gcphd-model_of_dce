library(mcce)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)
library(watercolours)
library(readr)
library(forcats)
library(tidyr)
library(pmwg)

args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 1) {
    stop("The model estimated data must be selected by name (old, std, reduced, reduced_more, perf1) as a passed argument", call.=FALSE)
}
model_run <- args[1]

# Load all the samples - only run one of the next four blocks
obj_name <- "sampler"
if (model_run == "std") {
  # Standard
  pref_file <- here::here("data", "output", "PrefDCE_m3vs9gjeqw6pmifp_OctEst.RData")
} else if (model_run == "reduced") {
  # Reduced
  pref_file <- here::here("data", "output", "PrefDCE_u2f0m2trer7kxcsj_OctEst.RData")
} else if (model_run == "reduced_more") {
  # Reduced More
  pref_file <- here::here("data", "output", "PrefDCE_j22mhqf8bl8ugaab_OctEst.RData")
} else if (model_run == "old") {
  # Older
  pref_file <- here::here("data", "output", "PrefDCE_2506730.rcgbcm_Estimation5Model.RData")
  obj_name <- "sampled"
} else if (model_run == "perf1") {
  # Standard
  pref_file <- here::here("data", "output", "EstimationPerformanceMarch", "PrefDCE_dvcljcisdh61vumc_TuffPerformance.RData")
} else {
  stop("Unknown estimated model selected")
}

pref_sampler <- get_samples(pref_file, final_obj = obj_name)
parameters <- pref_sampler$par_names
archs <- parameters[startsWith(parameters, "alpha")]
subjects <- pref_sampler$subjects
pref_samples <- extract_parameters(pref_sampler)

pdf(file = here::here("results", "Preferential", model_run, paste0("theta_mu_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
for (par in parameters) {
  g <- pref_samples %>%
    filter(subjectid == "theta_mu") %>%
    trace_plot(par) +
    scale_colour_watercolour()
  print(g)
  if (interactive()) {
    readline(prompt = "Press [enter] to continue")
  }
}
dev.off()

for (par in parameters) {
  print(paste("Generating pdf for", par))
  pdf(file = here::here("results", "Preferential", model_run, paste0(par, "_randeff_trace_", Sys.Date(), ".pdf")), width = 14.1, height = 7.53)
  for (subj_arr in split(subjects, ceiling(seq_along(subjects) / 4))) {
    g <- pref_samples %>%
      filter(subjectid %in% subj_arr) %>%
      trace_plot(par) +
      scale_colour_watercolour() +
      facet_wrap(~ subjectid, nrow = 2, ncol = 2)
    print(g)
  }
  dev.off()
}

model_medians <- pref_samples %>%
  filter(stageid == "sample", parameter %in% archs) %>%
  get_summary(tform = exp) %>%
  arch_medians()

model_medians %>%
  order_arch_medians() %>%
  arch_plot +
  scale_fill_watercolour()

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    model_run,
    paste0("EstimatedArch_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

par_medians <- pref_samples %>%
  filter(stageid == "sample", !(parameter %in% archs)) %>%
  get_summary()

par_colours <- c("t0" = "grey",
                 "A" = "grey", "b_acc" = "#73842E", "b_rej" = "#D0781C",
                 "v_acc_p_H" = "#569F72", "v_acc_p_L" = "#407755", "v_acc_p_D" = "#2B5039",
                 "v_rej_p_H" = "#F29F40", "v_rej_p_L" = "#EF8C1A", "v_rej_p_D" = "#BF6C0D",
                 "v_acc_r_H" = "#85C0FF", "v_acc_r_L" = "#47A0FF", "v_acc_r_D" = "#0A81FF",
                 "v_rej_r_H" = "#BBA9A0", "v_rej_r_L" = "#A2887C", "v_rej_r_D" = "#83685D")

par_medians %>%
  filter(subjectid != "theta_mu") %>%
  mutate(colour = par_colours[parameter]) %>%
  ggplot(aes(x = parameter, y = exp(value), fill = colour)) +
  geom_boxplot() +
  ylim(c(0, 7.5)) +
  scale_fill_identity()

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    model_run,
    paste0("EstimatedRE_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)

par_medians %>%
  filter(subjectid == "theta_mu") %>%
  mutate(colour = par_colours[parameter]) %>%
  ggplot(aes(x = parameter, y = exp(value), fill=colour)) +
  geom_col() +
  scale_fill_identity()

ggsave(
  filename = here::here(
    "results",
    "Preferential",
    model_run,
    paste0("EstimatedPars_", Sys.Date(), ".png")
  ),
  dpi = 200,
  width = 14.1,
  height = 7.53,
  units = "in"
)
