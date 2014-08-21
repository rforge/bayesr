library("BayesR")

data("marital.nz", package = "VGAM")
m <- model.matrix(~ -1 + mstatus, data = marital.nz)
colnames(m) <- gsub("/", "", colnames(m))
colnames(m) <- gsub("mstatus", "", colnames(m))
marital.nz <- cbind(marital.nz, m)


b0 <- bayesr(mstatus ~ age, family = multinomial, data = marital.nz,
  reference = "Married/Partnered", engine = "JAGS", n.iter = 1200, burnin = 200)

b1 <- bayesr(mstatus ~ sx(age), family = multinomial,
  data = marital.nz, reference = "Married/Partnered",
  engine = "BayesX", verbose = TRUE, n.iter = 1200, burnin = 200)

f <- list(
  DivorcedSeparated ~ sx(age),
  Single ~ sx(age),
  Widowed ~ sx(age)
)

b2 <- bayesr(f, family = multinomial, data = marital.nz,
  engine = "BayesX", verbose = FALSE, n.iter = 1200, burnin = 200)

f <- list(
  DivorcedSeparated ~ age,
  Single ~ age,
  Widowed ~ age
)

b3 <- bayesr(f, family = multinomial, data = marital.nz,
  engine = "JAGS", n.iter = 1200, burnin = 200)
