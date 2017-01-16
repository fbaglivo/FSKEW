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
    
    u_high=0.75; %Features
    b_high=0.3; 
    d_high=0.3;
    
    
    u_low=0.1; %Binding
    u_bkp(t)=u_low;
%     u_low=0.1;
    b_low=0.3;
    d_low=0.3;

    switch electrode_run_type
           
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
    
    
    Complete=[Bind' Feat'];
    
    Hfeat(t)=0.25*log(det(cov(zscore(Feat))));
    Hbind(t)=0.25*log(det(cov(zscore(Bind))));
    H(t)=0.5*log(det(cov(zscore(Complete'))));
    
end