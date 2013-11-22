library("BayesR")

data("rent", package = "gamlss.data")

rent$R2 <- rent$R / 1000

f <- R2 ~ Sp + Sm + B + H + L + loc + s(Fl) + s(A) 

b0 <- bayesr(f, family = gamma.BayesR, data = rent)

f <- R2 ~ sx(Fl) + sx(A)

b1 <- bayesx2(f, family = gamma.BayesR, data = rent)
