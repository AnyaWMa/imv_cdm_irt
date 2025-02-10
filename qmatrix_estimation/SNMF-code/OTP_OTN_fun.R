
OTP_OTN_fun<-function(trueQ,estQ){
  OTPtemp<-0
  OTNtemp<-0
  N<-nrow(trueQ)
  att<-ncol(trueQ)
  for (i in 1:N) {
    for (k in 1:att) {
        if(trueQ[i,k]==1 & estQ[i,k]==0){
          OTPtemp<-OTPtemp+1
        }
         if(trueQ[i,k]==0 & estQ[i,k]==1){
           OTNtemp<-OTNtemp+1
         } 
        
      }
  }
  oneSum<-sum(trueQ)
  zeroSum<-N*att-oneSum
  OTP<-OTPtemp/oneSum
  OTN<-OTNtemp/zeroSum
  result<-c(OTP,OTN)
  return(result)
  
}