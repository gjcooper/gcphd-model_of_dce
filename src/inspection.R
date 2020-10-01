library(coda)
require(tidyverse)
require(tcltk)
require(pmwg)

#Load in the data into the global environment
f <- tkgetOpenFile(
        title = "RData file",
        initialdir = here::here("data/output")) %>%
  as.character()

load(f)

# *Assumes* final object is named sampled
knitr::kable(table(sampled$samples$stage))

#' Plot the model parameter estimates from the sample stage.
#'
#' @param plot_obj The pmwgs object containing the estimates
#'
#' @return Nothing, side effect is it creates a plot
plot_theta_mu <- function(plot_obj) {
  samples <- as_mcmc(plot_obj, filter="sample")
  dimnames(samples) <- list(NULL,parameters)
  plot(samples, smooth=TRUE)
}

#' Plot the evidence for each model by subject
#'
#' @param plot_obj - the pmwg sampler object containing estimates
#' @param relative - Whether to plot the evidence as relative or absolute
#'
#' @return None - side effect is the creation of the plot
look_at_alphas <- function(plot_obj, relative=TRUE) {
  tmus <- as_mcmc(plot_obj, filter="sample")
  dimnames(tmus) <- list(NULL, parameters)
  tmu_alphas <- tmus[, 1:7] %>% data.frame() %>% tibble()
  tmu_alphas$subjectid <- 'theta_mu'

  res <- as_mcmc(plot_obj, selection="alpha", filter="sample") %>%
    lapply(FUN = function(X) {
      dimnames(X) <- list(NULL, parameters)
      data.frame(X[, 1:7])
    }) %>%
  bind_rows(.id="subjectid") %>%
  tibble() %>%
  rbind(tmu_alphas)


  medians <- res %>% group_by(subjectid) %>% summarise(across(.fns=median)) %>%
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

#Look at the estimates in absolute and relative terms, then use coda mcmc plot for theta_mu
look_at_alphas(sampled)
look_at_alphas(sampled, relative=FALSE)
plot_theta_mu(sampled)
