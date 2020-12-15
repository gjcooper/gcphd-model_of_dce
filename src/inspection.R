library(coda)
require(tidyverse)
require(tcltk)
require(pmwg)
require(mcmcplots)

get_data <- function(final_obj = 'sampled') {
  #Load in the data into the global environment
  f <- tkgetOpenFile(
          title = "RData file",
          initialdir = here::here("data/output")) %>%
    as.character()

  load(f, envir = (e <- new.env()))

  sampled <- e[[final_obj]]
  # *Assumes* final object is named sampled
  knitr::kable(table(sampled$samples$stage))
  sampled
}

#' Plot the model parameter estimates from the sample stage.
#'
#' @param plot_obj The pmwgs object containing the estimates
#'
#' @return Nothing, side effect is it creates a plot
plot_theta_mu <- function(plot_obj) {
  samples <- as_mcmc(plot_obj, filter="sample")
  dimnames(samples) <- list(NULL,plot_obj$par_names)
  plot(samples, smooth=TRUE)
}

collect_samples <- function(plot_obj) {
  tmus <- as_mcmc(plot_obj, filter="sample")
  dimnames(tmus) <- list(NULL, plot_obj$par_names)
  tmu_alphas <- tmus[, 1:7] %>% data.frame() %>% tibble()
  tmu_alphas$subjectid <- 'theta_mu'

  as_mcmc(plot_obj, selection="alpha", filter="sample") %>%
    lapply(FUN = function(X) {
      dimnames(X) <- list(NULL, plot_obj$par_names)
      data.frame(X[, 1:7])
    }) %>%
  bind_rows(.id="subjectid") %>%
  tibble() %>%
  rbind(tmu_alphas)
}

#' Plot the evidence for each model by subject
#'
#' @param plot_obj - the pmwg sampler object containing estimates
#' @param relative - Whether to plot the evidence as relative or absolute
#'
#' @return None - side effect is the creation of the plot
look_at_alphas <- function(plot_obj, relative=TRUE) {
  medians <- collect(plot_obj) %>%
    group_by(subjectid) %>%
    summarise(across(.fns=median)) %>%
    mutate_at(vars(contains("alpha")), exp) %>%
    pivot_longer(-subjectid) %>%
    rename(Parameter=name) %>%
    group_by(subjectid) %>%
    mutate(allvals = sum(value)) %>%
    mutate(relative_evidence = value/allvals) %>%
    mutate(Parameter = as.factor(Parameter))
  
  par_labels = sapply(strsplit(levels(medians$Parameter), "_"), function(x) {x[2]})

  if (relative) {
    ggplot(medians, aes(x=subjectid, y=relative_evidence, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels)
  } else {
    ggplot(medians, aes(x=subjectid, y=value, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels)
  }
}


get_medians <- function(res) {
  res %>%
    group_by(subjectid) %>%
    summarise(across(.fns = median)) %>%
    mutate_at(vars(contains("alpha")), exp) %>%
    pivot_longer(-subjectid) %>%
    rename(Parameter = name) %>%
    group_by(subjectid) %>%
    mutate(allvals = sum(value)) %>%
    mutate(relative_evidence = value / allvals) %>%
    mutate(Parameter = as.factor(Parameter))
}

compare <- function(original, recovery, relative = TRUE) {
  original <- collect_samples(original)
  recovery <- collect_samples(recovery)


  omedians <- get_medians(original)
  rmedians <- get_medians(recovery)
  par_labels = sapply(strsplit(levels(omedians$Parameter), "_"), function(x) {x[2]})

  medians <- bind_rows(estimates = omedians, recovery = rmedians, .id = "source")

  if (relative) {
    ggplot(medians, aes(x=subjectid, y=relative_evidence, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Relative Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels) +
      facet_grid(rows = vars(source))
  } else {
    ggplot(medians, aes(x=subjectid, y=value, fill=Parameter)) +
      geom_col() +
      xlab("Subject Identifier") +
      ylab("Absolute Evidence") +
      scale_fill_brewer(palette="Set2", name = "Model", labels = par_labels) +
      facet_grid(rows = vars(source))
  }
}

compare_data <- function(original, recovery, comparison = "rt") {
  odata <- original$data %>% as_tibble()
  rdata <- recovery$data %>% as_tibble()
  alldata <- bind_rows(original = odata, recovery = rdata, .id = "source")
  if (comparison == "rt") {
    alldata %>%
      filter(rt < 5) %>%
      ggplot(mapping = aes(x=rt, colour=source)) +
      geom_density() +
      facet_wrap(~ cell)
  } else if (comparison == "rtbox") {
    alldata %>%
      filter(rt < 5) %>%
      ggplot(mapping = aes(x=cell, y = rt, colour=source)) +
      geom_boxplot()
  } else if (comparison == "response") {
    alldata %>%
      ggplot(mapping = aes(x = response, group=source, fill=source)) +
      geom_bar(position="dodge") +
      facet_wrap(~ subject)
  }
}

#Look at the estimates in absolute and relative terms, then use coda mcmc plot for theta_mu
recovery <- get_data()
original <- get_data()
compare(original, recovery)
compare_data(original, recovery, "rt")

# Looking at Multichannel try3
load(here::here('data', 'output', 'Task1_MultiChannelTry3.RData'), envir = (try1_e <- new.env()))
sample_data <- try1_e$sampled
mcmcplot(as_mcmc(sample_data, filter="sample"))
mcmcplot(as_mcmc(sample_data, select='alpha', filter="sample"))
sample_df <- sample_data %>%
  as_mcmc() %>%
  as_tibble()

sample_df %>%
  select(starts_with('v_')) %>%
  pivot_longer(
    everything(),
    names_to = c('drift', 'response', 'attribute', 'salience'),
    names_transform = list(
      response = ~ readr::parse_factor(.x, levels=c('acc', 'rej')),
      attribute = ~ readr::parse_factor(.x, levels=c('p', 'r')),
      salience = ~ readr::parse_factor(.x, levels=c('H', 'L', 'D'))),
    names_sep="_") %>%
  select(-drift) %>%
  rename(drift=value) %>%
  mutate(
   response = fct_recode(response, Accept = 'acc', Reject = "rej"),
   attribute = fct_recode(attribute, Price = "p", Rating = "r")
  ) %>%
  ggplot(mapping=aes(x=salience, y=drift)) + geom_boxplot(aes(fill=salience)) + facet_grid(vars(response), vars(attribute))
