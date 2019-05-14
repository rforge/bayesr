## --- setup ---
set.seed(111)

## --- libraries ---
library("bamlss")

## --- data ---
load("ALDIS_ERA5_subset.rda")

## --- build formula ---
p <- names(d_train)
p <- p[p != "counts"]
p <- paste0("s(",p ,", bs = 'ps')")
p <- paste0(p, collapse = " + ")
f <- as.formula(paste("~", p))
f <- list(update(f, counts ~ .), f)

## --- stability selection ---
sel <- stabsel(f, data = d_train, family = "ztnbinom", q = 16, B = 100)

pdf("stabsel.pdf")
plot(sel, show = 25)
dev.off()

save(sel, file = "stabsel.rda")

## --- fit selected model ---
newf <- formula(sel)
fam <- family(sel)

b <- bamlss(newf, data = d_train,             # standard interface
  family = fam, binning = TRUE,               # general arguments
  optimizer = boost, maxit = 1000,            # boosting arguments
  thin = 5, burnin = 1000, n.iter = 6000      # sampler arguments
)

fit <- as.data.frame(predict(b, newdata = d_eval, type = "parameter"))
d_eval <- cbind(d_eval, fit)

save(b, d_eval, fam, sel, file = "full_model.rda")



