function signal=henongen_func(u,b,d,tipo)

% Programa para simular registros de EEG usando mapas de Henon
% ('Detection of coupling with linear and nonlinear synchronization 
%   measures for EEG'.Bakhshayesh, 2014)
% ('Learning driver-response relationships from synchronization
% patterns'.Quiroga, 2000)
%
% x(k+1)=1.4+bx(k-1)-x^2(k)
% y(k+1)=1.4+dy(k-1)-[ux(k)+(1-u)y(k)]y(k)
%
% u = 0 --> no coupling   u = 1 --> full coupling 
%
% b = 0.3 y d = 0.3 --> identical systems  
% b = 0.3 y d = 0.1 --> identical systems  
%
% Fabricio Baglivo 2016

x1=0;
x2=0;
% 
% u=0.6;  % Coupling Strengh
% 
% b=0.3;
% d=0.1;

x(1)=x1;
x(2)=x2;

y(1)=0;
y(2)=0;
y(3)=0;
y(4)=0;

% Signal calculation

for k=2:65000
   
    x(k+1)=1.4+b*x(k-1)-x(k)^2;
    y(k+1)=1.4+d*y(k-1)-[u*x(k)+(1-u)*y(k)]*y(k);
    
end

%Graphical verification
%%

% for k=50000:60000
% 
%     
%     scatter(x(k+1),x(k),'filled','r');
% %     scatter3(x(k),x(k+2),x(k+4),'filled','r');
%     hold on
%     
% end

%%

signal(:,1)=x(30000:65000);
signal(:,2)=y(30000:65000);

%%

save(['C:/Users/Fabri/Documents/Documentos/PhD/Papers/EnProceso/Metodos/Scripts/Mats/signals_' num2str(u) 'coupled' tipo '.mat'], 'signal','u','b','d');

end