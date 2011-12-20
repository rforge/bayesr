run.bayesx <- function(wait = TRUE, ...) 
{
  if (.Platform$OS.type == "windows") {
    cmd <- shQuote(system.file(package="bayesxsrc", "libs", .Platform$r_arch, "BayesX.exe"))
  } else {
    fn <- shQuote(system.file(package="bayesxsrc", "libs", .Platform$r_arch, "BayesX"))
    cmd <- paste(shQuote(file.path(R.home(),"bin","R")), "CMD", fn)
  }
  system(cmd, wait=wait, ...)
}
