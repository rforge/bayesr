library("BayesR")

data("marital.nz", package = "VGAM")

b0 <- bayesr(mstatus ~ s(age), family = multinomial,
  data = marital.nz, reference = "Married/Partnered",
  engine = "JAGS")

b1 <- bayesr(mstatus ~ sx(age), family = multinomial,
  data = marital.nz, reference = "Married/Partnered",
  engine = "BayesX")
