library("bamlss")

data("marital.nz", package = "VGAM")
m <- model.matrix(~ -1 + mstatus, data = marital.nz)
colnames(m) <- gsub("/", "", colnames(m))
colnames(m) <- gsub("mstatus", "", colnames(m))
marital.nz <- cbind(marital.nz, m)


b0 <- bamlss(mstatus ~ s(age), family = multinomial, data = marital.nz, reference = "Married/Partnered", sampler = NULL)

b1 <- bamlss(mstatus ~ sx(age), family = multinomial,
  data = marital.nz, reference = "Married/Partnered",
  engine = "BayesX", verbose = TRUE, n.iter = 1200, burnin = 200)


f <- list(
  DivorcedSeparated ~ sx(age),
  Single ~ sx(age),
  Widowed ~ sx(age)
)

b2 <- bamlss(f, family = multinomial.bamlss(link = "probit"), data = marital.nz,
  engine = "BayesX", verbose = FALSE, n.iter = 1200, burnin = 200)

f <- list(
  DivorcedSeparated ~ s(age),
  Single ~ s(age),
  Widowed ~ s(age)
)

b3 <- bamlss(f, family = multinomial, data = marital.nz,
  engine = "JAGS", n.iter = 1200, burnin = 200)

b4 <- bamlss(f, family = multinomial, data = marital.nz,
  method = c("backfitting", "MCMC"), update = "optim2",
  propose = "mvn", n.iter = 12000, burnin = 200, thin = 10)
