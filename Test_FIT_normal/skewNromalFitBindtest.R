# version 23112016


library("sn")
setwd("C:/Users/Fabri/Documents/Documentos/PhD/Papers/EnProceso/Metodos/FSKEW/")

file_fit<-paste("CSV/FitSize.csv",sep ="")
fit_size_csv = read.csv(file_fit,sep=",",header=FALSE)
fit_size=fit_size_csv[,1]

fit<-matrix(data = NA, nrow = 0, ncol = fit_size, byrow = TRUE, dimnames = NULL)

for(i in 1:length(dir("CSV/",pattern="Binding")))  
  {
    
  file<-paste("CSV/Binding.csv",sep ="")
  mydata = read.csv(file,sep=",",header=FALSE)
  formula<-paste("cbind(",paste(names(mydata),collapse = ","),")~1")
  m3=selm(formula,family="SN",data=mydata,method="MLE")
  fit<-rbind(fit,cbind(t(coef(m3,"DP")),m3@logL))
        
  }
write.table(fit, file = "skewNromalFitedDataBind.csv",row.names=FALSE,na="",col.names=FALSE, sep=",")

