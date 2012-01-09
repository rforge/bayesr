## BayesX REML testing
library("BayesXsrc")
reml <- run.bayesx("reml.prg", verbose = FALSE)
writeLines(reml$log)
