make.label <-
function(cx, xnam, dimx, vx)
{
  if(cx == "random.bayesx")
    label <- paste("r(", xnam, sep = "")
  if(cx == "sm.bayesx") {
    if(dimx > 1L)
      label <- paste("s(", xnam[1L], ",", xnam[2L], sep = "")
    else
      label <- paste("s(", xnam, sep = "")
  }
  if(cx == "mrf.bayesx")
    label <- paste("s(", xnam, sep = "")
  if(cx == "geo.bayesx")
    label <- paste("s(", xnam[1], sep = "")
  if(is.null(vx))
    label <- paste(label, ")", sep = "")
  else
    label <- paste(label, "):", vx, sep = "")

  return(label)
}

