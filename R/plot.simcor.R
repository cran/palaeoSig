plot.simcor <-
function(x,sim,...){
  hist(sim[,2], xlim=range(c(sim[,2],x$res$wc)), xlab="Weighted correlation", main="",col="grey80", border=NA,...)
  abline(v=x$res$wc,col=1)
  text(x$res$wc,par()$usr[4]*.9,label="pH", srt=90, pos=2)
  box()
}

