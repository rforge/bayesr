library("R2BayesX")

## installing the BayesX command-line binary from R
install.bayesx(inst.dir = "/path/to/bin", source.dir = NULL)

## important, before models can be fitted function bayesx()
## needs to know the location of the command-line binary of BayesX!
options(bayesx.bin = "/path/to/bin/BayesX")

## on Windows, after installion of the GUI version of BayesX,
## the command-line binary is in the folder 'commandline'
## and is named 'bayesx.exe', e.g. set
options(bayesx.bin = "C:/BayesX/commandline/bayesx.exe")

## code of the example
## load the data and bnd object
data("ZambiaNutrition", "ZambiaBnd", package = "R2BayesX")

## the model formula
f <- stunting ~ sx(agechild) + sx(mbmi) +
  sx(district, bs = "gk", map = ZambiaBnd, full = TRUE)

## estimate with BayesX from R
b <- bayesx(f, family = "gaussian", method = "MCMC",
  data = ZambiaNutrition)

## plot all terms
plot(b)

## with map
plot(b, map = ZambiaBnd)

## plot effect of mbmi only
plot(b, term = "sx(mbmi)")

## effect of agechild including partial residuals
plot(b, term = "sx(agechild)", residuals = TRUE,
  cex = 0.1, rug = FALSE)

## map effect plot of spatial effect with
## changing position of the legend
plot(b, term = "sx(district)", map = ZambiaBnd,
  swap = TRUE, pos = "right")

## view additional options for model terms
bayesx.term.options(bs = "ps", method = "MCMC")
bayesx.term.options(bs = "ps", method = "REML")
bayesx.term.options(bs = "ps", method = "STEP")

## illustration example
## set up the formula
f <- stunting ~ memployment + urban + gender + meducation +
  sx(mbmi) + sx(agechild) +
  sx(district, bs = "mrf", map = ZambiaBnd) + r(district)

## estimation
set.seed(321)
zm <- bayesx(f, family = "gaussian", method = "MCMC",
  data = ZambiaNutrition, iterations = 12000,
  burnin = 2000, step = 10)

## model summary
summary(zm)

## plot effects for mbmi and agechild
plot(zm, term = c("sx(mbmi)", "sx(agechild)"))

## plot kernel density estimates of the spatial
## and the random effect
plot(zm, term = c("sx(district)", "r(district)"))

## map effect plot of the spatial effect
plot(zm, term = "sx(district)", map = ZambiaBnd, swap = TRUE)

## map effect plot of the random effect with
## same color and legend scaling as for the
## structured spatial effect
range <- c(-0.32, 0.32)
plot(zm, term = "r(district)", map = ZambiaBnd,
  swap = TRUE, range = range, lrange = range)

## plot sampling paths of coefficients
## of term 'sx(mbmi)'
plot(zm, term = "sx(mbmi)", which = "coef-samples")

## plot autocorrellation function of term 'sx(mbmi)'
## and maximum autocorrellation of all parameters
## of the model
plot(zm, term = "sx(mbmi)", which = "var-samples", acf = TRUE)
plot(zm, which = "max-acf")

## extract samples for term 'sx(mbmi)'
## and plot with the coda package
samples.mbmi <- samples(zm, term = "sx(mbmi)")
library("coda")
plot(as.mcmc(samples.mbmi))
