## BayesX testing
testfiles <- c("mcmc.prg", "reml.prg", "step.prg", "BayesX-tests.R", "data.raw")

mcmc <- run.bayesx("mcmc.prg")
reml <- run.bayesx("reml.prg")
step <- run.bayesx("step.prg")

files <- list.files()
file.remove(files[!files %in% testfiles])
