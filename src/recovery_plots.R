library(mcce)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(pmwg)


#' Plot the evidence for each model by subject
#'
#' @param plot_obj - the pmwg sampler object containing estimates
#' @param relative - Whether to plot the evidence as relative or absolute
#'
#' @return None - side effect is the creation of the plot
plot_alphas <- function(plot_obj, relative = TRUE) {
  medians <- extract_alphas(plot_obj) %>%
    get_medians()

  if (relative) {
    ggplot(medians, aes(x = subjectid, y = rel_val, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model")
  } else {
    ggplot(medians, aes(x = subjectid, y = value, fill = Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette = "Set2", name = "Model")
  }
}



compare <- function(original, recovery, relative = TRUE) {
  original <- extract_alphas(original) %>%
    get_medians()
  recovery <- extract_alphas(recovery) %>%
    get_medians()

  medians <- bind_rows(
    estimates = original,
    recovery = recovery,
    .id = "source"
  )

}

compare_data <- function(original, recovery, comparison = "rt") {
  odata <- original$data %>%
    as_tibble()
  rdata <- recovery$data %>%
    as_tibble()
  alldata <- bind_rows(original = odata, recovery = rdata, .id = "source")
  if (comparison == "rt") {
    alldata %>%
      filter(rt < 5) %>%
      unite("cell", price:rating, sep='') %>%
      ggplot(mapping = aes(x = rt, colour = source)) +
      geom_density() +
      facet_wrap(~cell)
  } else if (comparison == "rtbox") {
    alldata %>%
      filter(rt < 5) %>%
      unite("cell", price:rating, sep='') %>%
      ggplot(mapping = aes(x = cell, y = rt, colour = source)) +
      geom_boxplot()
  } else if (comparison == "response") {
    alldata %>%
      mutate(accept=factor(accept, labels=c('reject', 'accept'))) %>%
      ggplot(mapping = aes(x = accept, group = source, fill = source)) +
      geom_bar(position = "dodge") +
      facet_wrap(~subject)
  }
}
# Look at the estimates in absolute and relative terms, then use coda mcmc plot for theta_mu


# Load all the samples
data_location <- here::here("data", "output", "5ModelRecovery")
recovery_files <- c(
  "NumericVDCE_CB_wvuhRfvyySjv_untagged.RData",
  "NumericVDCE_FPP_1878515.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_IEX_1878513.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_IST_1878369.rcgbcm_5ModelRecovery.RData",
  "NumericVDCE_MW_1878516.rcgbcm_5ModelRecovery.RData"
)

recovery_samples <- lapply(recovery_files, function(x) {
  get_samples(here::here(data_location, x))
})
names(recovery_samples) <- sapply(strsplit(recovery_files, "_"), "[[", 2)

original_samples <- get_samples(
  here::here(
    "data",
    "output",
    "NumericVDCE_1878182.rcgbcm_Estimation5Model.RData"
  )
)

all_samples <- c(recovery_samples, Original=original_samples)

recovery_medians <- sapply(recovery_samples, function(x) {
  extract_alphas(x) %>% get_medians()
  },
  USE.NAMES = TRUE,
  simplify = FALSE
)

original_medians <- original_samples %>%
  extract_alphas %>%
  get_medians

all_medians <- append(recovery_medians, original_medians)

model_order <- c("CB", "FPP", "MW", "IST", "IEX", "Original")

medians <- bind_rows(recovery_medians, .id = "source") %>%
  bind_rows(Original = original_medians) %>%
  replace_na(list(source = "Original")) %>%
  mutate(source = factor(source, levels = model_order))

ggplot(medians, aes(x = subjectid, y = rel_val, fill = Parameter)) +
  geom_col() +
  xlab("Subject Identifier") +
  ylab("Relative Evidence") +
  scale_fill_brewer(palette = "Dark2", name = "Model") +
  facet_grid(rows = vars(source))

dev.new()
ggplot(medians, aes(x = subjectid, y = value, fill = Parameter)) +
  geom_col() +
  xlab("Subject Identifier") +
  ylab("Absolute Evidence") +
  scale_fill_brewer(palette = "Dark2", name = "Model") +
  facet_grid(rows = vars(source))
