% This script takes binding a feature condition for faces stimulus during
% retrieval stage and calculate MI using Skewness technique. 
% The last 512 samples are selected for the analysis.
% Thi script works toghether with Skew fit R script
% Adapted by Fabricio Baglivo 2016 from Sergio Lew script.

clear all
close all
clc

stage='PostRetention'

load(['P12/' stage '.mat']);

BindinERPs_RED=cond(1).data(:,512:end,:);
FeaturesERPs_RED=cond(2).data(:,512:end,:);

nUnits=size(BindinERPs_RED,1);

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);
system('rm CSV/*');



%%
trials=size(BindinERPs_RED,3);

ok=0;

while ok==0

    vector=randperm(90);
    seleccion=vector(1:10);
    
    for iter=1:trials


        Bind(iter,:)=squeeze(mean(BindinERPs_RED(seleccion,:,iter),2));
        Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(seleccion,:,iter),2));


    end

    if (rank(Bind)==10) && (rank(Feat)==10);ok=1;end;
    
end

Complete=[Bind' Feat'];
csvwrite(['CSV/Complete.csv'],Complete');

idx=randperm(130);
Complete_Shuffle=Complete(:,idx);
Bind=Complete(:,1:65);
Feat=Complete(:,66:130);
csvwrite(['CSV/Binding.csv'],Bind');
csvwrite(['CSV/Features.csv'],Feat');



%%
%!unset DYLD_LIBRARY_PATH; Rscript  skewNromalFitBind.R

%%


!Rscript skewNromalFitBind_v17112016.R
snParam = csvread('skewNromalFitedDataBind.csv');
% snParam(snParam==-999999)=NaN;
!Rscript skewNromalFitFeat_v17112016.R
snParamBind = csvread('skewNromalFitedDataFeat.csv');
% snParamBind(snParam==-999999)=NaN;
!Rscript  skewNromalFit_v17112016.R
snParamFeat = csvread('skewNromalFitedData.csv');
% snParamFeat(snParamGO==-999999)=NaN;

%%

for n=1:size(snParam,1)
    
    alfa=snParam(n,66:75);
   
    init=11;
    for i=1:10 %for 10 channels
        
        l=10-i ;% #components    
        
        omega(i,i:10)=snParam(1,init:(init+l));
        omega(i:10,i)=snParam(1,init:(init+l));
    
        init=(init+l+1);
    
    end
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfa*alfa')*a>b))=a(find(sqrt(alfa*alfa')*a>b));
    W(find(sqrt(alfa*alfa')*a<=b))=-a(find(sqrt(alfa*alfa')*a<=b));
    H(n) = 1/2*log((det(omega))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
    
    %BIND
    
    alfaBind=snParamBind(n,66:75);
   
    init=11;
    for i=1:10 %for 10 channels
        
        l=10-i ;% #components    
        
        omegaBind(i,i:10)=snParamBind(1,init:(init+l));
        omegaBind(i:10,i)=snParamBind(1,init:(init+l));
    
        init=(init+l+1);
    
    end
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfaBind*alfaBind')*a>b))=a(find(sqrt(alfaBind*alfaBind')*a>b));
    W(find(sqrt(alfaBind*alfaBind')*a<=b))=-a(find(sqrt(alfaBind*alfaBind')*a<=b));
    HBind(n) = 1/2*log((det(omegaBind))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaBind*alfaBind')*W)));
    
    %FEAT
    
    alfaFeat=snParamFeat(n,66:75);
   
    init=11;
    for i=1:10 %for 10 channels
        
        l=10-i ;% #components    
        
        omegaFeat(i,i:10)=snParamFeat(1,init:(init+l));
        omegaFeat(i:10,i)=snParamFeat(1,init:(init+l));
    
        init=(init+l+1);
    
    end
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfaFeat*alfaFeat')*a>b))=a(find(sqrt(alfaFeat*alfaFeat')*a>b));
    W(find(sqrt(alfaFeat*alfaFeat')*a<=b))=-a(find(sqrt(alfaFeat*alfaFeat')*a<=b));
    HFeat(n) = 1/2*log((det(omegaFeat))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaFeat*alfaFeat')*W)));
        
    
end

MI=H - 1/2*HFeat - 1/2*HBind



