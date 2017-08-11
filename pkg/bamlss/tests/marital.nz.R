library("bamlss")

data("marital.nz", package = "VGAM")
m <- model.matrix(~ -1 + mstatus, data = marital.nz)
colnames(m) <- gsub("/", "", colnames(m))
colnames(m) <- gsub("mstatus", "", colnames(m))
marital.nz <- cbind(marital.nz, m)

b0 <- bamlss(mstatus ~ s(age), family = "multinomial",
  data = marital.nz, reference = "Married/Partnered", sampler = FALSE,
  optimizer = boostm, maxit = 1000)

f <- list(
  DivorcedSeparated ~ s(age),
  Single ~ s(age),
  Widowed ~ s(age)
)

b1 <- bamlss(f, family = "multinomial", data = marital.nz)

