% v24112016
%
%
% This script takes binding a feature condition for faces stimulus during
% retrieval stage and calculate MI using Skewness technique. 
% The last 512 samples are selected for the analysis.
% Thi script works toghether with Skew fit R script
%
% It is prepeared for n electrodes analysis using random electrodes
% or pre defined electrodes
% It can do Normal, trial shuffle or data shuffle run
%
% Adapted by Fabricio Baglivo 2016 from Sergio Lew script.

clear all
close all
clc

electrode_run_type='test'; %'random_electrodes'; %selected_electrodes 
data_run_type= 'data_shuffle'; %'normal';%trial_shuffle %data_shuffle

electrode_number=90;
electrode_selected_number=10; 

electrodes=[1 11 17 21 30 45 52];  %Case of selected electrodes

stage='PostRetention'

load(['P12/' stage '.mat']);

BindinERPs_RED=cond(1).data(:,512:end,:);
FeaturesERPs_RED=cond(2).data(:,512:end,:);

nUnits=size(BindinERPs_RED,1);

setenv('PATH', [getenv('PATH') ';C:\Program Files\R\R-3.3.2\bin\']);
system('rm CSV/*');



%%
trials=size(BindinERPs_RED,3);


for t=1
    
    u_high=0.4; %Features
    b_high=0.3; 
    d_high=0.3;
    
    
    u_low=0.1; %Binding
    u_bkp(t)=u_low;
%     u_low=0.1;
    b_low=0.3;
    d_low=0.3;

    switch electrode_run_type
        
        case 'random_electrodes'
            
            switch data_run_type
                
                case 'normal'
                    
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
                    
                    
                case 'trial_shuffle'
                    
                    ok=0;
                    
                    while ok==0
                        
                        vector=randperm(electrode_number);
                        selection=vector(1:electrode_selected_number);
                        
                        for iter=1:trials
                            
                            
                            Bind(iter,:)=squeeze(mean(BindinERPs_RED(selection,:,iter),2));
                            Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(selection,:,iter),2));
                            
                            
                        end
                        
                        if (rank(Bind)==10) && (rank(Feat)==10);ok=1;end;
                        
                    end
                    
                    
                    
                    %mix codition: Trial permutation
                    Complete=[Bind' Feat']';
                    idx=randperm(trials);
                    Complete_Shuffle=Complete(idx,:);
                    Bind=Complete(1:trials/2,:);
                    Feat=Complete(trials/2+1:trials,:);
                    
                    
                case 'data_shuffle'
                    
                    %data shuffle
                    
                    shufflevect=randperm(size(cond(1).data,2));
                    B=cond(1).data(:,shufflevect,:);
                    shufflevect=randperm(size(cond(2).data,2));
                    F=cond(1).data(:,shufflevect,:);
                    
                    BindinERPs_RED=B(:,512:end,:);
                    FeaturesERPs_RED=F(:,512:end,:);
                    
                    
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
                    
            end
            
        case 'selected_electrodes'
            
            
            electrode_selected_number=size(electrodes,2)
            selection=electrodes;
            
            switch data_run_type
                
                case 'normal'
                    
                    
                    
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
                    
                case 'trial_shuffle'
                    
                    
                    for iter=1:trials
                        
                        Bind(iter,:)=squeeze(mean(BindinERPs_RED(selection,:,iter),2));
                        Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(selection,:,iter),2));
                    end
                    
                    
                    
                    if (rank(Bind)==size(Bind,2)) && (rank(Feat)==size(Feat,2))
                        
                        ok=1;
                    else
                        
                        error('Non-independent data matrix');
                        
                    end
                    
                    
                    %mix codition: Trial permutation
                    Complete=[Bind' Feat']';
                    idx=randperm(trials);
                    Complete_Shuffle=Complete(idx,:);
                    Bind=Complete(1:trials/2,:);
                    Feat=Complete(trials/2+1:trials,:);
                    
                    
                case 'data_shuffle'
                    
                    %data shuffle
                    
                    shufflevect=randperm(size(cond(1).data,2));
                    B=cond(1).data(:,shufflevect,:);
                    shufflevect=randperm(size(cond(2).data,2));
                    F=cond(1).data(:,shufflevect,:);
                    
                    BindinERPs_RED=B(:,512:end,:);
                    FeaturesERPs_RED=F(:,512:end,:);
                    
                    
                    
                    for iter=1:trials
                        
                        
                        Bind(iter,:)=squeeze(mean(BindinERPs_RED(selection,:,iter),2));
                        Feat(iter,:)=squeeze(mean(FeaturesERPs_RED(selection,:,iter),2));
                        
                        
                    end
                    
                    
                    
                    if (rank(Bind)==size(Bind,2)) && (rank(Feat)==size(Feat,2))
                        
                        
                        ok=1;
                    else
                        
                        error('Non-independent data matrix');
                        
                    end
                    
            end
            
        case 'test'
            
            clear BindinERPs_RED FeaturesERPs_RED
            
            electrode_number=2;
            electrode_selected_number=2;
            
            %Bind
            
            tipo='low_conenction';
            u=u_low;
            b=b_low;d=d_low; % symetrical system
            signal=henongen_func(u,b,d,tipo);
            signal(:,1)=signal(:,1)-mean(signal(:,1));
            signal(:,2)=signal(:,2)-mean(signal(:,2));
             signal_low=signal;
            
            for i=1:1000%ceil(size(signal,1)/512)-1
                
                BindinERPs_RED(:,:,i)=signal((i-1)*32+1:i*32,:)';
                Bind(i,:)=squeeze(mean(BindinERPs_RED(:,:,i),2));
                
            end
            

%             Bind(:,1)=randn(65,1);
%             Bind(:,2)=randn(65,1);
              

     %FEAT
            
            tipo='hig_conenction';
            u=u_high;
            b=b_high;d=d_high; % Symetrical system
            signal=henongen_func(u,b,d,tipo);

            signal(:,1)=signal(:,1)-mean(signal(:,1));
            signal(:,2)=signal(:,2)-mean(signal(:,2));
            
            
            for i=1:1000%ceil(size(signal,1)/512)-1
                
                FeaturesERPs_RED(:,:,i)=signal((i-1)*32+1:i*32,:)';
%                 FeaturesERPs_RED(:,:,i)=signal(1:512,:)';
                Feat(i,:)=squeeze(mean(FeaturesERPs_RED(:,:,i),2));
                
            end
            
%             Feat(:,1)=rand(65,1);
%             Feat(:,2)=rand(65,1);

                            
    end
    
    % FIT vector:
    %
    %  n -> mean values
    %  sum([1:1:electrode_selected_number]) -> Omega Matrix
    %  n -> Aplha vector
    %  1 -> likelihood
    
    fit_size=2*electrode_selected_number+sum([1:1:electrode_selected_number])+1;
    
    %%
    
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
        
        likelihood(t)=snParam(end);
        
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
        
        
        likelihoodbind(t)=snParamBind(end);
        
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
        
        likelihoodFeat(t)=snParamFeat(end);
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

end

%%

save('Mat/Aasimetric.mat','HH','HHf','HHb');
    
%%

figure
subplot(1,2,1)
bar([H HBind HFeat],0.5)
set(gca,'xticklabels',{'H';'Hnc';'Hc'})
grid on
subplot(1,2,2)
bar(MI)
set(gca,'xticklabels',{'MI'})

figure;plot(Bind);xlabel('Low Connection')
figure;plot(Feat);xlabel('High Connection')

%%
figure
for i=1:9
    
    
    bar([HH(i) HHb(i) HHf(i)],0.5)
    set(gca,'xticklabels',{'H';'Hnc';'Hc'})
grid on

    pause(1)
    
end
