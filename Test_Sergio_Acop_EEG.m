clear all
close all
clc

electrode_run_type='test'; %'random_electrodes'; %selected_electrodes 
data_run_type= 'data_shuffle'; %'normal';%trial_shuffle %data_shuffle
% 
% electrode_number=90;
% electrode_selected_number=10; 

electrodes=[1 11];  %Case of selected electrodes

stage='Retention'

load(['P12/' stage '.mat']);

BindinERPs_RED_temp=cond(1).data(:,1:100,:); %512:end ; retention stage
FeaturesERPs_RED_temp=cond(2).data(:,1:100,:);

complete_temp=cat(3,BindinERPs_RED_temp,FeaturesERPs_RED_temp);

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);
system('rm CSV/*');



%%

BindinERPs_RED=BindinERPs_RED_temp;
FeaturesERPs_RED=FeaturesERPs_RED_temp;

for t=1
    
    %Shuffle
%     
%     shufflevect=randperm(size(complete_temp,3));
%     BindinERPs_RED=complete_temp(:,:,shufflevect(1:65));
%     FeaturesERPs_RED=complete_temp(:,:,shufflevect(65:130));
%     
    
    trials=size(BindinERPs_RED,3);
    electrode_selected_number=size(electrodes,2)
    selection=electrodes;
    nUnits=size(BindinERPs_RED,1);
    
    for iter=1:trials
        
        
        Bind(iter,:)=squeeze(mean(BindinERPs_RED(selection,:,iter),2));
        Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(selection,:,iter),2));
        
        
    end
    
    %Test independence bewteen signals
    
    if (rank(Bind)==size(Bind,2)) && (rank(Feat)==size(Feat,2))
        
        ok=1;
    else
        
        error('Non-independent data matrix');
        
    end
    
    
    Complete=[Bind' Feat'];
    
    Hfeatnorm=0.25*log(det(cov(zscore(Feat))));
    Hbindnorm=0.25*log(det(cov(zscore(Bind))));
    Hnorm=0.5*log(det(cov(zscore(Complete'))));
    
    
    fit_size=2*electrode_selected_number+sum([1:1:electrode_selected_number])+1;
    
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
    
    for n=1:size(snParam,1)
        
        alfainit=electrode_selected_number+sum([1:1:electrode_selected_number])+1;
        alfaend=alfainit+electrode_selected_number-1;
        
        alfa=snParam(n,alfainit:alfaend);
        
        likelihood=snParam(end);
        
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
        
        
        likelihoodbind=snParamBind(end);
        
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
        
        likelihoodFeat=snParamFeat(end);
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
    
    MI(t)=H - 1/2*HFeat - 1/2*HBind
    
    HH(t)=H;
    HHf(t)=HFeat;
    HHb(t)=HBind;


%%
end


Miprom=mean(MI);
HHprom=mean(HH);
HHfprom=mean(HHf);
HHbprom=mean(HHb);