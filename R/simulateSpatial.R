`simulateSpatial` <-
function(geo.xy, nsim=999, y,varmod, centre=F){
  if(any(duplicated(geo.xy)))stop("cannot have co-located samples")
  sims<-krige(sim~1,~x+y,data=geo.xy,dummy=T,nsim=nsim,beta=ifelse(centre,0,mean(y)),model=varmod,newdata=geo.xy)
  sims<-sims[,-(1:2)]
  if(centre) sims<-sweep(sims,2,mean)
  return(sims)
}

