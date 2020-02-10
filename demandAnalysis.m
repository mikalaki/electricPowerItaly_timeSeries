%% Ergasia Xronoseirwn : Zisou Charilaos AEM 9213 ,Karatzas Michalis AEM 9137
close all; clear; clc;

% import data
load('demandData.mat');

% %mikalaki code
% teamNumber=7;
% 
% %computing the time and the regionNumber ,we have to examine.
% time = mod(teamNumber,24) +1; %8
% regionNumber = mod(teamNumber,7) +1 +4; % column 5
% 
% %load the xls file data. 
% italyPowerData=xlsread('ElectricPowerItaly.xls','demand');
% 
% %getting the timeserie we want to examine 
% demand1=italyPowerData(italyPowerData(:,4)==8,5);
% 
% 
% %mikalaki end


%constants
per = 7;
maorder = 7;
maxtau = 100;
Tmax = 10;  
p=4;
q=1;
alpha = 0.05;
zalpha = norminv(1-alpha/2);
n = length(demand);

% plot unprocessed time series
figure(1)
clf
plot(demand)
%splot(demand,'-o')
hold on
xlabel('t')
ylabel('y(t)')
title('Unprocessed demand time series')

%% LINEAR 1
% % detrend time series
% maDemand = movingaveragesmooth2(demand, maorder);
% detrended = demand - maDemand;
% figure(2)
% plot(detrended)
% title(['Detrended demand time series by MA(', num2str(maorder), ') smoothing'])

% 1a.detrend time series First Differences
detrended=zeros(364,1);
for i=1:364
    detrended(i)=demand(i+1)-demand(i);
end
figure(2)
plot(detrended)
title(['Detrended demand time series by first differencies'])

% deseason time series
deseasonedc = seasonalcomponents(detrended, per);
%deseasonedma = movingaverageseasonal(detrended, per);
deseasoned =  detrended - deseasonedc;
figure(3)
plot(deseasoned)
title('Deseasoned demand time series')

%% LINEAR 2
% Autocorrelation of deseasoned time series using average seasonal component
figure(4)
ac1 = autocorrelation(deseasoned(~isnan(deseasoned)), maxtau);
autlim = zalpha/sqrt(n);
figure(5)
clf
hold on
for ii=1:maxtau
    plot(ac1(ii+1,1)*[1 1],[0 ac1(ii+1,2)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('r(\tau)')
title('Autocorrelation')

% Ljung-Box Portmanteau Test
figure(6)
tittxt = ('Deseasoned time series');
[h2V,p2V,Q2V] = portmanteauLB(ac1(2:maxtau+1,2),n,alpha,tittxt);

% Partial autocorrelation
pac1 = parautocor(deseasoned, maxtau);
figure(7)
clf
hold on
for ii=1:maxtau
    plot(ac1(ii+1,1)*[1 1],[0 pac1(ii)],'b','linewidth',1.5)
end
plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
xlabel('\tau')
ylabel('\phi_{\tau,\tau}')
title('Partial autocorrelation of deseasoned demand')

%% LINEAR 3
% %AIC criterion for different models( model with minimun AIC criterion is
% %the bestfit.
% akaike(deseasoned,Tmax);

% Fit ARMA Model
%[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(deseasoned, p, q, Tmax);
%aicS = akaike(deseasoned, Tmax);

%ARMA model parameters.
p=6
%p=3
q=0

%fitting the model
[nrmseV,phiV,thetaV,SDz,~,~,armamodel]=fitARMA(deseasoned, p, q, Tmax);
% 
%plotting the NRMSE for Tmax steps forward
hold off;
plot(nrmseV);
title(sprintf('NRMSE of ARMA(%d,%d) fit in demand time series, autocorrelation',p,q));

xlabel("steps forward");
ylabel("value");

%% LINEAR 4
%predictions=[];
indexes=[];
nrmse=[];
step=30;
for n1 = 20:step:n %6 different points in total
   indexes=[indexes n1];
   
   [nrmseV,~,~,~] = predictARMAnrmse(deseasoned,5,0,1,n-n1,[]);
   nrmse=[nrmse nrmseV];
   
   %[predictions] =[predictions predictARMAmultistep(deseasoned,n1,5,0,1,[])];
   %[preV] = predictARMAmultistep(xV,n1,p,q,Tmax,'example');
end

figure()

% plot(indexes,predictions,"-o");
% hold on;
plot(indexes,nrmse,"-o");
title("NRMSE for one step forward prediction with AR(5) model");
xlabel("Size of learning test");
legend("NRMSE");




