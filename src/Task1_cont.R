require(pmwg)
require(rtdists)
library(MCMCpack)

#Get input/output filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	exit("Input file must be entered to continue run")
}
outfile <- paste0("data/output/", args[1])

load(outfile)

adapted <- run_stage(burned, stage = "adapt", iter = 5000, particles = 500)

save.image(outfile)

sampled <- run_stage(adapted, stage = "sample", iter = 5000, particles = 100)

save.image(outfile)
