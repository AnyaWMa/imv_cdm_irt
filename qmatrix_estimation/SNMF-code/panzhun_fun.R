
panzhun_fun<-function(truealpha,estalpha){
  numP=0
  numPP=0
  N<-nrow(truealpha)
  att<-ncol(truealpha)
  for (i in 1:N) {
    t=0
    for (k in 1:att) {
      if (estalpha[i,k]==truealpha[i,k]){
        numP=numP+1
        t=t+1
      }
    }
    if (t==att)
      numPP=numPP+1
  }
  MMR=numP/(N*att)
  PMR=numPP/N
  result<-c(MMR,PMR)
  return(result)
  
}