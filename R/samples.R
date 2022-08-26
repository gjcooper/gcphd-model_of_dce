tmu <- function(sampler, func = median) {
  sampler %>%
    as_mcmc(filter = "sample") %>%
    as_tibble() %>%
    setNames(sampler$par_names) %>%
    summarise_all(func)
}

alph <- function(sampler, func = median) {
  sampler %>%
    as_mcmc(selection = "alpha", filter = "sample") %>%
    sapply(function(x) {
             x %>%
               as_tibble() %>%
               setNames(sampler$par_names) %>%
               summarise_all(func)},
           simplify=FALSE) %>%
    bind_rows(.id = "SubjectID")
}
