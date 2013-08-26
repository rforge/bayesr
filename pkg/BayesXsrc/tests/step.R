## BayesX STEP testing
library("BayesXsrc")
step <- run.bayesx("step.prg", verbose = FALSE)
fx1 <- read.table("step_f_x1_pspline.res", header = TRUE)
fx4 <- read.table("step_f_x4_pspline.res", header = TRUE)
print(fx1)
print(fx4)
