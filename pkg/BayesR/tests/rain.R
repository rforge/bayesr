library("BayesR")

if(file.exists("~/tmp/homstart.rda")) {
  load("~/tmp/homstart.rda")
} else {
  dpath <- system.file(package = "BayesR", "data")
  if(!file.exists(file <- file.path(dpath, "homstart.rda"))) {
    homstart_data(dir = dirname(file), load = TRUE)
  } else load(file)
}

rain <- subset(homstart, year >= 2008 & raw > 0)

f <- list(
  sqrt(raw) ~ s(day, bs = "cc") + s(elevation),
  ~ s(day, bs = "cc") + s(elevation)
)

b1 <- bayesr(f, data = rain, method = "MP2", family = gF(truncreg, point = 0), n.samples = 50)


library("truncreg")

set.seed(1071)
n <- 10000
sigma <- 4
alpha <- 2
beta <- 1
x <- rnorm(n, mean = 0, sd = 2)
eps <- rnorm(n, sd = sigma)
y <- alpha + beta * x + eps
d <- data.frame(y = y, x = x)

## truncated response
d$yt <- ifelse(d$y > 1, d$y, NA)

b0 <- truncreg(yt ~ x, data = d, point = 1, direction = "left")

b <- bayesr(yt ~ x, data = d, family = gF(truncreg, point = 1))
mf <- model.frame(b)
eta <- fitted(b)
dy <- gF(truncreg, point = 1)$d(mf$yt, eta)
hist(mf$yt, freq = FALSE)
lines(dy[order(mf$yt)] ~ mf$yt[order(mf$yt)], col = 2)

