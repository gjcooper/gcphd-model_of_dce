sampler <- pmwgs(
  data = model_data,
  pars = parameters,
  ll_func = dirichlet_func,
  prior = priors
)

cores <- vars$cores

sampler <- init(sampler, particles = 10)

sampler <- run_stage(sampler, stage = "burn", iter = 5, particles = 10, n_cores = cores)

load(testfile, envir = e <- new.env())

print("Testing equality of new run to saved tests")
new_theta_mu <- sampler$samples$theta_mu
old_theta_mu <- e$sampler$samples$theta_mu

if (identical(old_theta_mu, new_theta_mu)) {
  print("Old and new are identical")
} else {
  print("Old and new are not identical")
  all.equal(old_theta_mu, new_theta_mu)
  save.image(outfile)
}
