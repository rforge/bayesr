findBayesX <- function()
{
  if(.Platform$OS.type == "windows") {
    path <- shQuote(system.file(package="BayesXsrc", "libs", .Platform$r_arch, "BayesX.exe"))
  } else {
    path <- shQuote(system.file(package="BayesXsrc", "libs", .Platform$r_arch, "BayesX"))
  }
  if(path == "") {
    path <- getOption("bayesx.bin")
  }
  return(path)
}
