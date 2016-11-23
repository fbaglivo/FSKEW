library("sn")
setwd("C:/Users/Fabri/Documents/Documentos/PhD/Papers/EnProceso/Metodos/FSKEW/")
fit<-matrix(data = NA, nrow = 0, ncol = 8, byrow = TRUE, dimnames = NULL)

for(i in 1:length(dir("CSV/",pattern="FR1FR2_pair_Bind_")))  
  {
    
  file<-paste("CSV/FR1FR2_pair_Bind_",i,".csv",sep ="")
  mydata = read.csv(file,sep=",",header=FALSE,col.names=1:2)
  m3=selm(cbind(X1,X2)~1,family="SN",data=mydata,method="MLE")
  fit<-rbind(fit,cbind(t(coef(m3,"DP")),m3@logL))
        
  }
write.table(fit, file = "skewNromalFitedDataBind.csv",row.names=FALSE,na="",col.names=FALSE, sep=",")

