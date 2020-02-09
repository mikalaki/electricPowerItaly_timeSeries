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
hold on
xlabel('t')
ylabel('y(t)')
title('Unprocessed demand time series')

% detrend time series
maDemand = movingaveragesmooth2(demand, maorder);
detrended = demand - maDemand;
figure(2)
plot(detrended)
title(['Detrended demand time series by MA(', num2str(maorder), ') smoothing'])

% deseason time series
deseasonedc = seasonalcomponents(detrended, per);
%deseasonedma = movingaverageseasonal(detrended, per);
deseasoned =  detrended - deseasonedc;
figure(3)
plot(deseasoned)
title('Deseasoned demand time series')

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
title('Partial autocorrelation')


% Fit ARMA Model
%[nrmseV,phiallV,thetaallV,SDz,aicS,fpeS]=fitARMA(deseasoned, p, q, Tmax);
%aicS = akaike(deseasoned, Tmax);



