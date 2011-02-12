obscor.sig <-
function(x,sim){
  res<-(sum(x$res$wc<=sim[,2])+1)/(nrow(sim)+1)
  return(res)
}

