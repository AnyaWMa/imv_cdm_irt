#####packages and functions used 
library(GDINA)
library(CDM)
library(NMF)
source('utility.R')
source('compute_rss_fun.R')

####Fractional subtraction data
data(data.fraction4)
R<-data.fraction4$data

J<-ncol(R)
N<-nrow(R)
JJ<-(J+1)/2+1

rss<-compute_rss_fun(R,JJ)  

#####calculate the st based on rss########################
  st<-rep(0,JJ)
  f<-rss
  ep<-0.5     #当分母为0时的因子
  for (j in 2:(JJ-1)) {
      if((f[j]-f[j+1])==0){
        st[j]<-abs(f[j]-f[j-1])/ep
      }else{
        st[j]<-abs(f[j]-f[j-1])/abs(f[j+1]-f[j])
      }
  }
  print(which.max(st))

######保存结果##########################
write.csv(st,file = 'snmf_K3_st.csv',row.names = FALSE)
