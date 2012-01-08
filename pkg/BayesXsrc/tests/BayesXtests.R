## BayesX testing
BayesXtests <- function()
{
  prg <- c("mcmc.prg", "reml.prg", "step.prg", "BayesX-tests.R", "data.raw")
  if (.Platform$OS.type == "windows") {
    bin <- system.file(package = "BayesXsrc", "libs", .Platform$r_arch, "BayesX.exe")
  } else {
    bin <- system.file(package = "BayesXsrc", "libs", .Platform$r_arch, "BayesX")
  }

  ## mcmc
  mcmc <- system(paste(bin, prg[1L]))
  if(mcmc == 0) {
    files <- list.files()
    files <- files[!files %in% prg]
    file.remove(files)
  } else stop("problems testing method MCMC!")

  ## reml
  reml <- system(paste(bin, prg[2L]))
  if(reml == 0) {
    files <- list.files()
    files <- files[!files %in% prg]
    file.remove(files)
  } else stop("problems testing method REML!")

  ## step
  step <- system(paste(bin, prg[3L]))
  if(step == 0) {
    files <- list.files()
    files <- files[!files %in% prg]
    file.remove(files)
  } else stop("problems testing method STEP!")

  cat("everything is fine : )!\n")
}
