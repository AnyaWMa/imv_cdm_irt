library(CDM)
library(NMF)
source('qmatrix_estimation/SNMF-code/utility.R')
source('qmatrix_estimation/SNMF-code/compute_rss_fun.R')

dataset <- redivis::user("datapages")$dataset("item_response_warehouse",version='v2.0')

dataset_tables <- dataset$list_tables()

df <- dataset$table("content_literacy_intervention")$to_data_frame()
ntimes<-2

x0 <- df %>% 
  arrange(item)
resp0<-data.frame(irw::long2resp(x0))
id<-resp0$id
resp0$id<-NULL

resp0

# determine K
R <- resp0
J<-ncol(R)
N<-nrow(R)
JJ<-(J+1)/2+1

rss<-compute_rss_fun(R,JJ)  

#####calculate the st based on rss########################
st<-rep(0,JJ)
f<-rss
ep<-0.5     #????ĸΪ0ʱ??????
for (j in 2:(JJ-1)) {
  if((f[j]-f[j+1])==0){
    st[j]<-abs(f[j]-f[j-1])/ep
  }else{
    st[j]<-abs(f[j]-f[j-1])/abs(f[j+1]-f[j])
  }
}
print(which.max(st))


# q_matrix identification
K<-2    
R0<- resp0
R<-del.zeros(R0)  
print(nrow(R))

######estimate the Q-matrix used the SNMF method
time1=as.POSIXlt(Sys.time()) 
R.nmf<- nmf(t(!R),K,method = 'snmf/l') 
time2=as.POSIXlt(Sys.time())
use.time=difftime(time2,time1,units="secs")
print(use.time)     

R.W<-R.nmf@fit@W

R.W

R.W_converted <- ifelse(R.W > 1, 1, 0)

R.W_converted