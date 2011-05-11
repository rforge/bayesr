mum <-
function(x)
{
  rn <- rownames(x)
  a <- duplicated(x[,1L], fromLast = FALSE)
  b <- duplicated(x[,1L], fromLast = TRUE)
  a[a != b] <- TRUE
  if(any(tot <- grepl(":total", rn))) {
    x <- x[!tot,]
    rn <- rn[!tot]
  }
  if(any(a)) {
    drn <- rn[a]
    nc <- nchar(drn)
    rmrn <- drn[nc < max(nc, na.rm = TRUE)]
    x <- x[!(rn %in% rmrn),]
  }
    
  return(x)
}

