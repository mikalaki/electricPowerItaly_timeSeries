%% Ergasia Xronoseirwn : Zisou Charilaos AEM 9213 ,Karatzas Michalis AEM 9137
close all; clear; clc;

% import data
load('pricesData.mat');

% %mikalaki code
% teamNumber=7;
% 
% %computing the time and the regionNumber ,we have to examine.
% time = mod(teamNumber,24) +1; %8
% regionNumber = mod(teamNumber,7) +1 +4; % column 5
% 
% %load the xls file data. 
% italyPowerData=xlsread('ElectricPowerItaly.xls','prices');
% 
% %getting the timeserie we want to examine 
% prices=italyPowerData(italyPowerData(:,4)==8,5);
% %mikalaki end


%constants
per = 7;
maorder = 7;
maxtau = 100;
alpha = 0.05;
Tmax=20;
zalpha = norminv(1-alpha/2);
n = length(prices);

% plot unprocessed time series
figure(1)
clf
% plot(prices,'-o')
plot(prices)
hold on
xlabel('t')
ylabel('Prices')
title('Unprocessed prices time series')

%% LINEAR 1 
% % detrend time series MA
% maPrices = movingaveragesmooth2(prices, maorder);
% detrended = prices - maPrices;
% figure(2)
% plot(detrended)
% title(['Detrended prices time series by MA(', num2str(maorder), ') smoothing'])

% 1a.detrend time series First Differences
detrended=zeros(364,1);
for i=1:364
    detrended(i)=prices(i+1)-prices(i);
end
figure(2)
plot(detrended)
title(['Detrended prices time series by first differencies'])

%1b. deseason time series
deseasonedc = seasonalcomponents(detrended, per);
%deseasonedma = movingaverageseasonal(detrended, per);
deseasoned =  detrended - deseasonedc;
figure(3)
plot(deseasoned)
title('Deseasoned (and detrended) prices time series')


%%  LINEAR 2
%2. Autocorrelation of deseasoned time series using average seasonal component
ac1 = autocorrelation(deseasoned, maxtau);
autlim = zalpha/sqrt(n);
figure(4)
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
title('Autocorrelation of deseasoned time series')
% tittxt = sprintf('deseasoned time series by ave seasonal comp (%d), Ljung-Box test',per);
% figure(5)

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

%% LINEAR 3
% %getting the AICS for different models, in order to choose
% akaike(deseasoned,Tmax);
% %result is p=3 and q=3

%ARMA model parameters.
p=3
%p=3
q=3

%fitting the model
[nrmseV,phiV,thetaV,SDz,aicS,fpeS,armamodel]=fitARMA(deseasoned, p, q, Tmax);

%printing the coefficients
fprintf('Estimated coefficients of phi(B):\n');
disp(phiV')
fprintf('Estimated coefficients of theta(B):\n');
disp(thetaV')
fprintf('SD of noise: %f \n',SDz);

%Some criteria values 
fprintf('AIC: %f \n',aicS);
fprintf('FPE: %f \n',fpeS);

%plotting the NRMSE for Tmax steps forward
figure();
plot(nrmseV,"-o");
title("NRMSE");
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
title("NRMSE for one step forward prediction with AR(5) model, in prices time series");
xlabel("Size of learning test");
legend("NRMSE");


