## BayesX MCMC testing
library("BayesXsrc")
mcmc <- run.bayesx("mcmc.prg", verbose = FALSE)
fx1 <- read.table("mcmc_f_x1_pspline.res", header = TRUE)
fx2 <- read.table("mcmc_f_x2_pspline.res", header = TRUE)
print(head(fx1), digits = 4)
print(head(fx2), digits = 4)
