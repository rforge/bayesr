rmf <-
function(x) 
{
  for(i in 1L:length(x))
    for(char in c("+", "-", "*", ":", "^", "/", " ")) 
      x[i] <- gsub(char, "_", x[i], fixed = TRUE)

  return(x)
}

