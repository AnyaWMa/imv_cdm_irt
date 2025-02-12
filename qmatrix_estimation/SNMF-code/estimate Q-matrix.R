#####packages and functions used 
library(GDINA)
library(CDM)
library(NMF)
source('panzhun_fun.R') #
source('OTP_OTN_fun.R')   #
source('utility.R')

#############
 ########Fractional subtraction data£¬which is used in Xu&Zhang,2018;Chen et al.,2015
data(data.fraction4)
R0<-data.fraction4$data
 ####### the Q-matrix used in Xu&Zhang,2018
true_Q3p=c(1,0,0,
           1,0,0,
           1,0,0,
           0,1,0,
           0,0,1,
           0,0,1,
           0,0,1,
           1,0,1,
           1,1,0,
           1,1,1,
           1,1,0,
           0,1,0,
           0,1,0,
           0,1,1,
           1,1,1,
           0,1,0,
           1,1,1)
 true_Q3p=matrix(true_Q3p,3,17)
 q.matrix=t(true_Q3p)
###
 K<-ncol(q.matrix)     
 R<-del.zeros(R0)  
 print(nrow(R))

 ######estimate the Q-matrix used the SNMF method
 time1=as.POSIXlt(Sys.time()) 
 R.nmf<- nmf(t(!R),K,method = 'snmf/l') 
 time2=as.POSIXlt(Sys.time())
 use.time=difftime(time2,time1,units="secs")
 print(use.time)     
 
 R.W<-R.nmf@fit@W
 Q.d<-discreQ_fun(q.matrix,R.W)  
 
 ###
 Q_panzhun<-panzhun_fun(q.matrix,Q.d)
 EMR.Q<-Q_panzhun[1]
 PMR.Q<-Q_panzhun[2]
 
 OTP_OTN_Q<-OTP_OTN_fun(q.matrix,Q.d)
 OTP.Q<-OTP_OTN_Q[1]
 OTN.Q<-OTP_OTN_Q[2]
 
 output<-c(EMR.Q,PMR.Q,OTP.Q,OTN.Q)
 print(output)
 write.csv(Q.d,file = "SNMF-Q.csv")

 