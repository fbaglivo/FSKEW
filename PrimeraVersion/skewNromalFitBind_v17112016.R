library("sn")
setwd("C:/Users/Fabri/Documents/Documentos/PhD/Papers/EnProceso/Metodos/FSKEW/PrimeraVersion/")
fit<-matrix(data = NA, nrow = 0, ncol = 76, byrow = TRUE, dimnames = NULL)

for(i in 1:length(dir("CSV/",pattern="Binding")))  
  {
    
  file<-paste("CSV/Binding.csv",sep ="")
  mydata = read.csv(file,sep=",",header=FALSE)
  formula<-paste("cbind(",paste(names(mydata),collapse = ","),")~1")
  m3=selm(formula,family="SN",data=mydata,method="MLE")
  fit<-rbind(fit,cbind(t(coef(m3,"DP")),m3@logL))
        
  }
write.table(fit, file = "skewNromalFitedDataBind.csv",row.names=FALSE,na="",col.names=FALSE, sep=",")

