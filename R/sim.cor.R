sim.cor <-
function(spp,env,fos, n=99, ord=rda){
  res<-replicate(n,{
    mod<-WA(spp,runif(length(env)))
    pred<-predict(mod, fos)$fit[,1]
    RDA<-ord(fos~pred)
    optima<-mod$coef
    sco<-scores(RDA, display="spec", choice=1)
    abun<-t(t(colMeans(spp)))
    
    optima<-optima[intersect(rownames(optima),rownames(sco)),, drop=FALSE]
    abun<-abun[intersect(rownames(abun),rownames(sco)),, drop=FALSE]
    sco<-sco[intersect(rownames(sco),rownames(optima)),, drop=FALSE]
    x<-data.frame(optima, sco, abun=abun)    
    c(cr=abs(cor(x)[1,2]),wcr=abs(cov.wt(x[,1:2], wt=sqrt(x$abun), cor=TRUE)$cor[1,2]))
  
  })
  return(t(res))
}

