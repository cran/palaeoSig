plot.palaeoSig <-
function(x, names, pos=2,...){
  if(missing(names))names=names(x$EX)
  with(x,{
    hist(sim.ex, breaks=seq(min(sim.ex),max(sim.ex),length=20),xlim=c(0,MAX[1]*1.1), main="", xlab="Proportion variance explained", col="grey80", border=NA,...)
    abline(v=MAX, col=1, lwd=2, lty=3)
    abline(v=EX, col=1)     
    text(EX,par()$usr[4]*.9,label=names, srt=90, pos=pos)
  })
}

