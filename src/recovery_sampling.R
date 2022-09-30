sampler <- pmwgs(
  data = test_data,
  pars = parameters,
  ll_func = dirichlet_func,
  prior = priors
)

cores <- vars$cores

sampler <- init(sampler)

sampler <- run_stage(sampler, stage = "burn", iter = 5000, particles = 500, n_cores = cores)

save.image(outfile)

sampler <- run_stage(sampler, stage = "adapt", iter = 10000, particles = 500, n_cores = cores, n_unique = 40)

save.image(outfile)

sampler <- run_stage(sampler, stage = "sample", iter = 10000, particles = 100, n_cores = cores, pdist_update_n = NA)

save.image(outfile)
