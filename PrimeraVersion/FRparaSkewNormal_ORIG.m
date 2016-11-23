clear all
% close all
load data_good_4000_7400ms.mat

addpath('C:\Program Files\R\R-3.2.2\bin\')

w=320;
for iter=1:1
    idxt=1;
%    for t=[3500 3820 4170 4310]-w/2
     for t=(3500:50:5500)-w/2
        idxp=1;
        system('rm FRskew/*');
        for s=1:30
            nUnits=length(dataPFC.session(s).unit);
            tgoC=data_beh.session(s).trialGOc;
            tnogoC=data_beh.session(s).trialNOGOc;
            tgoI=data_beh.session(s).trialGOi;
            tnogoI=data_beh.session(s).trialNOGOi;
            tnogo=[tnogoC];
            tgo=[tgoC];
            
            tonos=[];
            tonos(tgo,1)=1;
            tonos(tnogo,1)=0;
            
            indTrials=sort([tgo,tnogo]);
            
            pGo=length(tgo)/length(indTrials);
            pNoGo=length(tnogo)/length(indTrials);
            MIaux=[];
            for u1=1:nUnits
                raster_basal=full(dataPFC.session(s).unit(u1).all(indTrials,3000:4000));
                raster1=full(dataPFC.session(s).unit(u1).all(indTrials,t:t+w));
                FR1=sum(raster1,2);%-mean(sum(raster_basal(:,1:5),2));
                raster1go=full(dataPFC.session(s).unit(u1).all(tgo,t:t+w));
                FR1go=sum(raster1go,2);%-mean(sum(raster_basal(:,1:5),2));
                raster1nogo=full(dataPFC.session(s).unit(u1).all(tnogo,t:t+w));
                FR1nogo=sum(raster1nogo,2);%-mean(sum(raster_basal(:,1:5),2));
                for u2=u1+1:nUnits
                    raster_basal=full(dataPFC.session(s).unit(u2).all(indTrials,3000:4000));
                    raster2=full(dataPFC.session(s).unit(u2).all(indTrials,t:t+w));
                    FR2=sum(raster2,2);%-mean(sum(raster_basal(:,1:5),2));
                    FR=[FR1 FR2];
                    raster2go=full(dataPFC.session(s).unit(u2).all(tgo,t:t+w));
                    FR2go=sum(raster2go,2);%-mean(sum(raster_basal(:,1:5),2));
                    FRgo=[FR1go FR2go];
                    raster2nogo=full(dataPFC.session(s).unit(u2).all(tnogo,t:t+w));
                    FR2nogo=sum(raster2nogo,2);%-mean(sum(raster_basal(:,1:5),2));
                    FRnogo=[FR1nogo FR2nogo];
                    csvwrite(['FRskew/FR1FR2_pair_' num2str(idxp) '.csv'],FR);
                    csvwrite(['FRskew/FR1FR2_pair_GO_' num2str(idxp) '.csv'],FRgo);
                    csvwrite(['FRskew/FR1FR2_pair_NOGO_' num2str(idxp) '.csv'],FRnogo);
                    idxp=idxp+1;
                end
            end
        end
        
        !unset DYLD_LIBRARY_PATH; Rscript  skewNromalFit.R
        snParam = csvread('skewNromalFitedData.csv');
        snParam(snParam==-999999)=NaN;
        !unset DYLD_LIBRARY_PATH; Rscript  skewNromalFitGO.R
        snParamGO = csvread('skewNromalFitedDataGO.csv');
        snParamGO(snParamGO==-999999)=NaN;
        !unset DYLD_LIBRARY_PATH; Rscript skewNromalFitNOGO.R
        snParamNOGO = csvread('skewNromalFitedDataNOGO.csv');
        snParamNOGO(snParamNOGO==-999999)=NaN;
        
        for n=1:size(snParam,1)
            alfa=snParam(n,6:7);
            omega=[snParam(n,3) snParam(n,4);snParam(n,4) snParam(n,5)];
                        
           
            M=100000;
            a=randn(100000,1);
            b=randn(100000,1);
            W(find(sqrt(alfa*alfa')*a>b))=a(find(sqrt(alfa*alfa')*a>b));
            W(find(sqrt(alfa*alfa')*a<=b))=-a(find(sqrt(alfa*alfa')*a<=b));
            H(n) = 1/2*log((det(omega))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
%             H(n) =  - mean(2*log(normcdf(sqrt(alfa*alfa')*W)));
%             H(n) = 1/2*log((det(omega)));

            alfaGO=snParamGO(n,6:7);
            omegaGO=[snParamGO(n,3) snParamGO(n,4);snParamGO(n,4) snParamGO(n,5)];
            M=100000;
            a=randn(100000,1);
            b=randn(100000,1);
            W(find(sqrt(alfaGO*alfaGO')*a>b))=a(find(sqrt(alfaGO*alfaGO')*a>b));
            W(find(sqrt(alfaGO*alfaGO')*a<=b))=-a(find(sqrt(alfaGO*alfaGO')*a<=b));
             HGO(n) = 1/2*log((det(omegaGO))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaGO*alfaGO')*W)));
%             HGO(n) = - mean(2*log(normcdf(sqrt(alfaGO*alfaGO')*W)));
%             HGO(n) = 1/2*log((det(omegaGO)));
        
            alfaNOGO=snParamNOGO(n,6:7);
            omegaNOGO=[snParamNOGO(n,3) snParamNOGO(n,4);snParamNOGO(n,4) snParamNOGO(n,5)];
            M=100000;
            a=randn(100000,1);
            b=randn(100000,1);
            W(find(sqrt(alfaNOGO*alfaNOGO')*a>b))=a(find(sqrt(alfaNOGO*alfaNOGO')*a>b));
            W(find(sqrt(alfaNOGO*alfaNOGO')*a<=b))=-a(find(sqrt(alfaNOGO*alfaNOGO')*a<=b));
            HNOGO(n) = 1/2*log((det(omegaNOGO))) + 1 + log(2*pi) - mean(2*log(normcdf(sqrt(alfaNOGO*alfaNOGO')*W)));
%             HNOGO(n) = - mean(2*log(normcdf(sqrt(alfaNOGO*alfaNOGO')*W)));
%             HNOGO(n) = 1/2*log((det(omegaNOGO)));
        end
        MI(idxt,:)=H - pGo * HGO - pNoGo * HNOGO;
        idxt=idxt+1
    end
end
figure;errorbar((3500:50:5500)-w/2,nanmean(MI'),nanstd(MI')/sqrt(size(MI',1)),'k')