sampler <- pmwgs(
  data = model_data %>% ungroup() %>% data.frame(),
  pars = parameters,
  ll_func = dirichlet_func,
  prior = priors
)

cores <- vars$cores

if (!exists("start_points")) {
  start_points <- NULL
}
sampler <- init(sampler, start_mu = start_points)

sampler <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = cores)

save.image(outfile)

sampler <- run_stage(sampler, stage = "adapt", iter = 10000, particles = 500, n_cores = cores)

save.image(outfile)

sampler <- run_stage(sampler, stage = "sample", iter = 10000, particles = 100, n_cores = cores)

save.image(outfile)




# Old ver: 14.5s
# Removed grouped_df: 7s
# Removed tibbles: 4.8s
# Removed do.call: 6.7s 5.6s 5.9s 8s 8.1s 5.7s 4.9s 7.7s 6.3s
#* Replace do.call: 5.4s 5.1s 6.1s 6.5s 6.7s 4.8s .3s 5.1s 5.3s Selected this
