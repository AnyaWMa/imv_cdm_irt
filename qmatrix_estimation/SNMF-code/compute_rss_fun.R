compute_rss_fun<-function(R,knum){
  #knum<-JJ
  Rt<-del.zeros(R)
  res<-nmf(t(Rt),seq(1,knum),method = 'snmf/l',.opt='vp5')  #
  onerss<-res$measures$rss
  return(onerss)
}

