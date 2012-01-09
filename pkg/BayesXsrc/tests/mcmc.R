## BayesX MCMC testing
library("BayesXsrc")
mcmc <- run.bayesx("mcmc.prg", verbose = FALSE)
writeLines(mcmc$log)
