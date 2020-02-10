%% Ergasia Xronoseirwn : Zisou Charilaos AEM 9213, Karatzas Michalis AEM 9137
%% Linear Analysis
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
maxtau = 100;
Tmax = 10;  
p=4;
q=1;
alpha = 0.05;
zalpha = norminv(1-alpha/2);
n = length(demand)-1;

% plot unprocessed time series
figure(1)
clf
plot(demand)
hold on
xlabel('t')
ylabel('y(t)')
title('Unprocessed demand time series')

%% LINEAR 1
% % Detrend time series by first order differences
detrended=zeros(364,1);
for i=1:364
    detrended(i)=demand(i+1)-demand(i);
end
figure(2)
plot(detrended)
title('Detrended demand time series by first differencies')

% deseason time series
deseasonedc = seasonalcomponents(detrended, per);
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
hold off

% Partial autocorrelation of deseasoned time series
pac1 = parautocor(deseasoned, maxtau);
figure(6)
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
title('Partial autocorrelation of demand time series')
hold off

% Ljung-Box Portmanteau Test
figure(7)
tittxt = ('Deseasoned time series');
[h2V,p2V,Q2V] = portmanteauLB(ac1(2:maxtau+1,2),n,alpha,tittxt);


%% LINEAR 3

% AIC
%aicMatrix = akaike(deseasoned, Tmax);
% BEST FIT IS ARMA(6,0)
p=6;
q=0;

% Fit ARMA Model
[nrmseARMA,~,~,~,~,~,armamodel]=fitARMA(deseasoned, p, q, Tmax);

figure(9)
plot(nrmseARMA,'-o')
ylabel("NRMSE");
xlabel("T");
title('Demand: NRMSE of ARMA(3,3) for T=1,2...10')

%% LINEAR 4
index=1;
nrmseAR5 = zeros(5,1);
for i=70:5:90
    [nrmseAR5(index),~,~,~] = predictARMAnrmse(deseasoned, 5, 0, 1, round(0.01*(100-i)*n), '');
    index=index+1;
end
figure(10)
plot((0.7:0.05:0.9).',nrmseAR5,'-o')
ylabel("NRMSE");
xlabel("Percent of training data");
title('Demand: NRMSE of AR(5) using 70%, 75%, 80%, 85%, 90% of training data')



