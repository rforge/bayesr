library("R2BayesX")


## Installing from R
install.bayesx(inst.dir = "/path/to/bin", source.dir = NULL)

## Important, before models can be fitted function bayesx()
## needs to know the location of the command-line binary of BayesX!
options(bayesx.bin = "/path/to/bin/BayesX")

## On Windows after installion of the GUI version of BayesX, the
## command-line binary is in the folder 'commandline'
## and is named 'bayesx.exe', e.g. set
options(bayesx.bin = "C:/BayesX/commandline/bayesx.exe")


## Code of the example
## Load the data and bnd object
data("ZambiaNutrition", "ZambiaBnd", package = "R2BayesX")

## The model formula
f <- stunting ~ sx(agechild) + sx(mbmi) +
  sx(district, bs = "gk", map = ZambiaBnd, full = TRUE)

## Estimate with BayesX from R
b <- bayesx(f, family = "gaussian", method = "MCMC",
  data = ZambiaNutrition)

## Plot all terms
plot(b)

## With map
plot(b, map = ZambiaBnd)

## Plot effect of mbmi only
plot(b, term = "sx(mbmi)")

## Effect of agechild including partial residuals
plot(b, term = "sx(agechild)", residuals = TRUE,
  cex = 0.1, rug = FALSE)

## Map effect plot of spatial effect with
## changing position of the legend
plot(b, term = "sx(district)", map = ZambiaBnd,
  swap = TRUE, pos = "right")


## View additional options for model terms
bayesx.term.options(bs = "ps", method = "MCMC")


## Illustration example
## Set up the formula
f <- stunting ~ memployment + urban + gender + meducation +
  sx(mbmi) + sx(agechild) +
  sx(district, bs = "mrf", map = ZambiaBnd) + r(district)

## Estimation
zm <- bayesx(f, family = "gaussian", method = "MCMC",
  iterations = 12000, burnin = 2000, step = 10,
  seed = 123, data = ZambiaNutrition)

## Model summary
summary(zm)


## Plot effects for mbmi and agechild
plot(zm, term = c("sx(mbmi)", "sx(agechild)"))

## Plot kernel density estimates of the spatial
## and the random effect
plot(zm, term = c("sx(district)", "r(district)"))

## Map effect plot of the spatial effect
plot(zm, term = "sx(district)", map = ZambiaBnd, swap = TRUE)

## Map effect plot of the random effect with
## same color and legend scaling as for the
## structured spatial effect
range <- c(-0.32, 0.32)
plot(zm, term = "r(district)", map = ZambiaBnd,
  swap = TRUE, range = range, lrange = range)

## Plot sampling paths of coefficients
## of term 'sx(mbmi)'
plot(zm, term = "sx(mbmi)", which = "coef-samples")

## Plot autocorrellation function of term 'sx(mbmi)'
## and maximum autocorrellation of all parameters
## of the model
plot(zm, term = "sx(mbmi)", which = "var-samples", acf = TRUE)
plot(zm, which = "max-acf")

