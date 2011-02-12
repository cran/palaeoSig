plot.obscor <-
function(x, xlab="WA optima", ylab="RDA scores", f=1,...){
  plot(x=x$x$Optima,y=x$x$RDA1, cex=x$x$abun*f, xlab=xlab, ylab=ylab,...)
}

