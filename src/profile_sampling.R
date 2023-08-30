Rprof()

sampler <- pmwgs(
  data = model_data,
  pars = parameters,
  ll_func = dirichlet_func,
  prior = priors
)

cores <- vars$cores

if (!exists("start_points")) {
  start_points <- NULL
}
sampler <- init(sampler, start_mu = start_points, particles = 10)

sampler <- run_stage(sampler, stage = "burn", iter = 20, particles = 100, n_cores = cores)

save.image(outfile)

sampler <- run_stage(sampler, stage = "adapt", iter = 100, particles = 100, n_cores = cores)

save.image(outfile)

Rprof(NULL)

summaryRprof()
