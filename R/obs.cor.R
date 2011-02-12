obs.cor <-
function(spp,env,fos, ord=rda){
  mod<-WA(spp,env)
  pred<-predict(mod, fos)$fit[,1]
  RDA<-ord(fos~pred)
  optima<-mod$coef
  sco<-scores(RDA, display="spec", choice=1)
  abun<-t(t(colMeans(spp)))
  
  optima<-optima[intersect(rownames(optima),rownames(sco)),, drop=FALSE]
  abun<-abun[intersect(rownames(abun),rownames(sco)),, drop=FALSE]
  sco<-sco[intersect(rownames(sco),rownames(optima)),, drop=FALSE]
  x<-data.frame(optima, sco, abun=abun)
  res<-list(x=x,res=list(wc=abs(cov.wt(x[,1:2], wt=sqrt(x$abun), cor=TRUE)$cor[1,2]),cc=abs(cor(x[,1:2])[1,2])))
  class(res)<-"obscor"
  return(res)
}

