x2int <- f2int <- function(x) 
{
  return(as.integer(as.numeric(as.character(x))))
}


# f2int <- function(x) 
# {
#   nl <- 1L:nlevels(x)
#   levels(x) <- nl
#
#   return(as.integer(paste(x)))
# }
