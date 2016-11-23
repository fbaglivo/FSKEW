% v23112016
%
%
% This script takes binding a feature condition for faces stimulus during
% retrieval stage and calculate MI using Skewness technique. 
% The last 512 samples are selected for the analysis.
% Thi script works toghether with Skew fit R script
% Adapted by Fabricio Baglivo 2016 from Sergio Lew script.

clear all
close all
clc

electrode_run_type='random_electrodes'; %selected_electrodes 
data_run_type='normal';%trial_shuffle %data_shuffle

stage='PostRetention'

load(['P12/' stage '.mat']);

BindinERPs_RED=cond(1).data(:,512:end,:);
FeaturesERPs_RED=cond(2).data(:,512:end,:);

nUnits=size(BindinERPs_RED,1);

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);
system('rm CSV/*');


%%

electrode_number=90;
electrode_selected_number=10;

% FIT vector:
%
%  n -> mean values
%  sum([1:1:electrode_selected_number]) -> Omega Matrix 
%  n -> Aplha vector
%  1 -> likelihood 

fit_size=2*electrode_selected_number+sum([1:1:electrode_selected_number])+1;

%%
trials=size(BindinERPs_RED,3);

ok=0;

while ok==0

    vector=randperm(electrode_number);
    selection=vector(1:electrode_selected_number);
    
    for iter=1:trials


        Bind(iter,:)=squeeze(mean(BindinERPs_RED(selection,:,iter),2));
        Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(selection,:,iter),2));


    end

    %Test independence bewteen signals
    if (rank(Bind)==electrode_selected_number) && (rank(Feat)==electrode_selected_number);ok=1;end;
    
end

csvwrite(['CSV/Binding.csv'],Bind);
csvwrite(['CSV/Features.csv'],Feat);
Complete=[Bind' Feat'];
csvwrite(['CSV/Complete.csv'],Complete');

csvwrite(['CSV/FitSize.csv'],fit_size);

%%
%!unset DYLD_LIBRARY_PATH; Rscript  skewNromalFitBind.R

%%


!Rscript skewNromalFitBind.R
snParam = csvread('skewNromalFitedDataBind.csv');
% snParam(snParam==-999999)=NaN;
!Rscript skewNromalFitFeat.R
snParamBind = csvread('skewNromalFitedDataFeat.csv');
% snParamBind(snParam==-999999)=NaN;
!Rscript  skewNromalFit.R
snParamFeat = csvread('skewNromalFitedData.csv');
% snParamFeat(snParamGO==-999999)=NaN;

%%

for n=1:size(snParam,1)
    
    alfainit=electrode_selected_number+sum([1:1:electrode_selected_number])+1;
    alfaend=alfainit+electrode_selected_number-1;
    
    alfa=snParam(n,alfainit:alfaend);
   
    init=electrode_selected_number+1;
    for i=1:electrode_selected_number %for 10 channels
        
        l=electrode_selected_number-i ;% #components    
        
        omega(i,i:electrode_selected_number)=snParam(1,init:(init+l));
        omega(i:electrode_selected_number,i)=snParam(1,init:(init+l));
    
        init=(init+l+1);
    
    end
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfa*alfa')*a>b))=a(find(sqrt(alfa*alfa')*a>b));
    W(find(sqrt(alfa*alfa')*a<=b))=-a(find(sqrt(alfa*alfa')*a<=b));
    H(n) = 1/2*log((det(omega))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
    
    %BIND
    
    alfaBind=snParamBind(n,alfainit:alfaend);
   
    init=electrode_selected_number+1;
    for i=1:electrode_selected_number 
        
        l=electrode_selected_number-i ;% #components    
        
        omegaBind(i,i:electrode_selected_number)=snParamBind(1,init:(init+l));
        omegaBind(i:electrode_selected_number,i)=snParamBind(1,init:(init+l));
    
        init=(init+l+1);
    
    end
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfaBind*alfaBind')*a>b))=a(find(sqrt(alfaBind*alfaBind')*a>b));
    W(find(sqrt(alfaBind*alfaBind')*a<=b))=-a(find(sqrt(alfaBind*alfaBind')*a<=b));
    HBind(n) = 1/2*log((det(omegaBind))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaBind*alfaBind')*W)));
    
    %FEAT
    
    alfaFeat=snParamFeat(n,alfainit:alfaend);
   
    init=electrode_selected_number+1;
    for i=1:electrode_selected_number %for 10 channels
        
        l=electrode_selected_number-i ;% #components    
        
        omegaFeat(i,i:electrode_selected_number)=snParamFeat(1,init:(init+l));
        omegaFeat(i:electrode_selected_number,i)=snParamFeat(1,init:(init+l));
    
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



