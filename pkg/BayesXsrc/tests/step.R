## BayesX STEP testing
library("BayesXsrc")
step <- run.bayesx("step.prg", verbose = FALSE)
writeLines(step$log)
