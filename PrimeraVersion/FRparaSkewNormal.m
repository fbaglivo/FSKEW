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

trial=1;

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);
system('rm CSV/*');

for iter=1:1
    idtxt=1;
    %    for t=[3500 3820 4170 4310]-w/2
    for u1=1:nUnits
        for u2=u1+1:nUnits
            
            FRbind=[squeeze(BindinERPs_RED(u1,:,trial))' squeeze(BindinERPs_RED(u2,:,trial))'];
            size(FRbind);
            FRfeat=[squeeze(FeaturesERPs_RED(u1,:,trial))' squeeze(FeaturesERPs_RED(u2,:,trial))'];
            csvwrite(['CSV/FR1FR2_pair_Feat_' num2str(idtxt) '.csv'],FRbind);
            csvwrite(['CSV/FR1FR2_pair_Bind_'  num2str(idtxt) '.csv'],FRfeat);
            idtxt=1+idtxt;
        end
    end
end

%!unset DYLD_LIBRARY_PATH; Rscript  skewNromalFitBind.R

%%

%!Rscript skewNromalFitBind.R
snParam = csvread('skewNromalFitedDataBind.csv');
% snParam(snParam==-999999)=NaN;
%!Rscript  skewNromalFitFeat.R
snParamGO = csvread('skewNromalFitedDataFeat.csv');
% snParamGO(snParamGO==-999999)=NaN;

%%

for n=1:size(snParam,1)
    alfa=snParam(n,6:7);
    omega=[snParam(n,3) snParam(n,4);snParam(n,4) snParam(n,5)];
    
    
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfa*alfa')*a>b))=a(find(sqrt(alfa*alfa')*a>b));
    W(find(sqrt(alfa*alfa')*a<=b))=-a(find(sqrt(alfa*alfa')*a<=b));
    %             H(n) = 1/2*log((det(omega))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
    %             H(n) =  - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
    H(n) = 1/2*log((det(omega)));
    
    alfaGO=snParamGO(n,6:7);
    omegaGO=[snParamGO(n,3) snParamGO(n,4);snParamGO(n,4) snParamGO(n,5)];
    M=100000;
    a=randn(100000,1);
    b=randn(100000,1);
    W(find(sqrt(alfaGO*alfaGO')*a>b))=a(find(sqrt(alfaGO*alfaGO')*a>b));
    W(find(sqrt(alfaGO*alfaGO')*a<=b))=-a(find(sqrt(alfaGO*alfaGO')*a<=b));
    %             HGO(n) = 1/2*log((det(omegaGO))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaGO*alfaGO')*W)));
    %             HGO(n) = - mean(2*log(normcdf(sqrt(alfaGO*alfaGO')*W)));
    HGO(n) = 1/2*log((det(omegaGO)));
    
%     alfaNOGO=snParamNOGO(n,6:7);
%     omegaNOGO=[snParamNOGO(n,3) snParamNOGO(n,4);snParamNOGO(n,4) snParamNOGO(n,5)];
%     M=100000;
%     a=randn(100000,1);
%     b=randn(100000,1);
%     W(find(sqrt(alfaNOGO*alfaNOGO')*a>b))=a(find(sqrt(alfaNOGO*alfaNOGO')*a>b));
%     W(find(sqrt(alfaNOGO*alfaNOGO')*a<=b))=-a(find(sqrt(alfaNOGO*alfaNOGO')*a<=b));
%     %             HNOGO(n) = 1/2*log((det(omegaNOGO))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaNOGO*alfaNOGO')*W)));
%     %             HNOGO(n) = - mean(2*log(normcdf(sqrt(alfaNOGO*alfaNOGO')*W)));
%     HNOGO(n) = 1/2*log((det(omegaNOGO)));
end
MI(idtxt,:)=H - HGO;
idtxt=idtxt+1


figure;errorbar((3500:50:5500)-w/2,nanmean(MI'),nanstd(MI')/sqrt(size(MI',1)),'k')