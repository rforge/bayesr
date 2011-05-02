library("R2BayesX")

## load forest health data and tree location map
data("ForestHealth")
data("BeechBnd")

## estimate model without spatial effect
fm1 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  s(age, bs = "ps") + s(inclination, bs = "ps") + 
  s(canopy, bs = "ps") + s(year, bs = "ps") + 
  s(elevation, bs = "ps"),
  family = "cumlogit", method = "REML", data = ForestHealth)
summary(fm1)
plot(fm1, term = c("s(age)", "s(inclination)", 
  "s(canopy)", "s(year)", "s(elevation)"))

## now include spatial effect
## warning: long runtime
fm2 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  s(age, bs = "ps") + s(inclination, bs = "ps") + 
  s(canopy, bs = "ps") + s(year, bs = "ps") + 
  s(elevation, bs = "ps") + 
  s(id, bs = "gk", xt = list(map = BeechBnd, full = TRUE)),
  family = "cumlogit", method = "REML", data = ForestHealth)
summary(fm2)
plot(fm2, term = c("s(age)", "s(inclination)", 
  "s(canopy)", "s(year)", "s(elevation)", "s(id)"),
  map = BeechBnd, pos = "topleft")

## compare effects for elevation and inclination
plot(c(fm1, fm2), term = c("s(elevation)", "s(inclination)"))
