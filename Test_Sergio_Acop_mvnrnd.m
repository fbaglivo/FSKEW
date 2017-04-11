clear all
close all
clc

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);

%%

electrode_selected_number=[1:1:10];
%%

indice=1;
for t=electrode_selected_number
    
    
    
    Mu=[zeros(10000,t)];
    SIGMA=diag([ones(t,1)],0);
    Rb=mvnrnd(Mu,SIGMA);
    
    Mu=[ones(10000,t)];
    Mu=Mu*5;
    SIGMA=diag([ones(t,1)],0);
%     
%     SIGMA(1,2)=0;
%     SIGMA(2,1)=0;
    
    Rf=mvnrnd(Mu,SIGMA);
    
    
    Bind=Rb;
    Feat=Rf;
        
    Complete=[Bind' Feat'];
    
    fit_size=2*t+sum([1:1:t])+1;
    
    csvwrite(['CSV/Binding.csv'],Bind);
    csvwrite(['CSV/Features.csv'],Feat);
    Complete=[Bind' Feat'];
    csvwrite(['CSV/Complete.csv'],Complete');
    
    csvwrite(['CSV/FitSize.csv'],fit_size);
    
    %%
    %!unset DYLD_LIBRARY_PATH; Rscript  skewNromalFitBind.R
    
    !Rscript skewNromalFitBind.R
    snParamBind = csvread('skewNromalFitedDataBind.csv');
    % snParam(snParam==-999999)=NaN;
    !Rscript skewNromalFitFeat.R
    snParamFeat = csvread('skewNromalFitedDataFeat.csv');
    % snParamBind(snParam==-999999)=NaN;
    !Rscript  skewNromalFit.R
    snParam = csvread('skewNromalFitedData.csv');
    % snParamFeat(snParamGO==-999999)=NaN;
    
    %%

    size(snParam,1)
    
    for n=1:size(snParam,1)
        
        alfainit=t+sum([1:1:t])+1;
        alfaend=alfainit+t-1;
        
        alfa=snParam(n,alfainit:alfaend);
        
        likelihood=snParam(end);
        
        init=t+1;
        for i=1:t %for 10 channels
            
            l=t-i ;% #components
            
            omega(i,i:t)=snParam(1,init:(init+l));
            omega(i:t,i)=snParam(1,init:(init+l));
            
            init=(init+l+1);
            
        end
        
        M=100000;
        a=randn(100000,1);
        b=randn(100000,1);
        W(find(sqrt(alfa*alfa')*a>b))=a(find(sqrt(alfa*alfa')*a>b));
        W(find(sqrt(alfa*alfa')*a<=b))=-a(find(sqrt(alfa*alfa')*a<=b));
        H(n) = 1/2*log((det(omega))) +  t/2*(1 + log(2*pi)) - mean(log(2*normcdf(sqrt(alfa*alfa')*W)));
        
        NoNorm(t)=H(n);
        
        sknormBind(t)= mean(log(2*normcdf(sqrt(alfa*alfa')*W)));
        
        %BIND
        
        alfaBind=snParamBind(n,alfainit:alfaend);
        
        
        likelihoodbind=snParamBind(end);
        
        init=t+1;
        for i=1:t
            
            l=t-i ;% #components
            
            omegaBind(i,i:t)=snParamBind(1,init:(init+l));
            omegaBind(i:t,i)=snParamBind(1,init:(init+l));
            
            init=(init+l+1);
            
        end
        
        M=100000;
        a=randn(100000,1);
        b=randn(100000,1);
        W(find(sqrt(alfaBind*alfaBind')*a>b))=a(find(sqrt(alfaBind*alfaBind')*a>b));
        W(find(sqrt(alfaBind*alfaBind')*a<=b))=-a(find(sqrt(alfaBind*alfaBind')*a<=b));
        HBind(n) = 1/2*log((det(omegaBind))) +  t/2*(1 + log(2*pi)) - mean(log(2*normcdf(sqrt(alfaBind*alfaBind')*W)));
        
        BindNoNorm(indice)=HBind(n);
        
        sknormBind(indice)= mean(log(2*normcdf(sqrt(alfaBind*alfaBind')*W)));
        
        %FEAT
        
        alfaFeat=snParamFeat(n,alfainit:alfaend);
        
        likelihoodFeat=snParamFeat(end);
        init=t+1;
        for i=1:t %for 10 channels
            
            l=t-i ;% #components
            
            omegaFeat(i,i:t)=snParamFeat(1,init:(init+l));
            omegaFeat(i:t,i)=snParamFeat(1,init:(init+l));
            
            init=(init+l+1);
            
        end
        
        M=100000;
        a=randn(100000,1);
        b=randn(100000,1);
        W(find(sqrt(alfaFeat*alfaFeat')*a>b))=a(find(sqrt(alfaFeat*alfaFeat')*a>b));
        W(find(sqrt(alfaFeat*alfaFeat')*a<=b))=-a(find(sqrt(alfaFeat*alfaFeat')*a<=b));
        HFeat(n) = 1/2*log((det(omegaFeat))) + t/2*(1 + log(2*pi)) - mean(log(2*normcdf(sqrt(alfaFeat*alfaFeat')*W)));
        
        %%%%% FALTA EL K/2!!!! ( 1 + log(2*pi))*K/2
        
        FeatNoNorm(indice)=HFeat(n);
        sknormFeat(indice)= mean(log(2*normcdf(sqrt(alfaFeat*alfaFeat')*W)));
    end
    
    MI(indice)=H - 1/2*HFeat - 1/2*HBind;
    
    HH(indice)=H;
    HHf(indice)=HFeat;
    HHb(indice)=HBind;

    indice=indice+1;

end

%%
% 
% Miprom=mean(MI);
% HHprom=mean(HH);
% HHfprom=mean(HHf);
% HHbprom=mean(HHb);