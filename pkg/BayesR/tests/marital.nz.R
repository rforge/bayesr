library("BayesR")

data("marital.nz", package = "VGAM")

b0 <- bayesr(mstatus ~ s(age), family = multinomial.BayesR(link = "probit"),
  data = marital.nz, reference = "Married/Partnered")

b1 <- bayesx2(mstatus ~ sx(age), family = multinomial.BayesR(link = "probit"),
  data = marital.nz, reference = "Married/Partnered")
