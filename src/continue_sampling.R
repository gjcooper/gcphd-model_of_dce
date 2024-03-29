cores <- vars$cores

if ("burn" %in% stages_to_run) {
  sampler <- run_stage(sampler, stage = "burn", iter = 2000, particles = 500, n_cores = cores)
}

save.image(outfile)

if ("adapt" %in% stages_to_run) {
  sampler <- run_stage(sampler, stage = "adapt", iter = 10000, particles = 500, n_cores = cores)
}

save.image(outfile)

if ("sample" %in% stages_to_run) {
  sampler <- run_stage(sampler, stage = "sample", iter = 10000, particles = 100, n_cores = cores)
}

save.image(outfile)
