library("sn")
setwd("/Users/sergiolew/Documents/sergio/Implante/DiscriminacionHeadFixedING/paper2014/scriptsMI/")
fit<-matrix(data = NA, nrow = 0, ncol = 8, byrow = TRUE, dimnames = NULL)

for(i in 1:length(dir("FRskew/",pattern="NOGO")))  
{
  file<-paste("FRskew/FR1FR2_pair_",i,".csv",sep ="")
  mydata = read.csv(file,sep=",",col.names=1:2)
  if(mean(mydata[,"X1"])>0.01 & mean(mydata[,"X2"])>0.01)
  {
    m3=selm(cbind(X1,X2)~1,family="SN",data=mydata,method="MLE")
    fit<-rbind(fit,cbind(t(coef(m3,"DP")),m3@logL))
  }
  else
  {
    fit<-rbind(fit,cbind(-999999,-999999,-999999,-999999,-999999,-999999,-999999,-999999))
  }
}
write.table(fit, file = "skewNromalFite)