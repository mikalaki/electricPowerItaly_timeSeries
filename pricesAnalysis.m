%% Ergasia Xronoseirwn : Zisou Charilaos AEM 9213 ,Karatzas Michalis AEM 9137
close all; clear; clc;

% % import data
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
% prices1=italyPowerData(italyPowerData(:,4)==8,5);
% 
% 
% %mikalaki end

% 
% %constants
% per = 7;
% maorder = 7;
% maxtau = 100;
% alpha = 0.05;
% zalpha = norminv(1-alpha/2);
% n = length(prices);
% 
% % plot unprocessed time series
% figure(1)
% clf
% plot(prices)
% hold on
% xlabel('t')
% ylabel('y(t)')
% title('Unprocessed prices time series')
% 
% % detrend time series
% maPrices = movingaveragesmooth2(prices, maorder);
% detrended = prices - maPrices;
% figure(2)
% plot(detrended)
% title(['Detrended prices time series by MA(', num2str(maorder), ') smoothing'])
% 
% % deseason time series
% deseasonedc = seasonalcomponents(detrended, per);
% %deseasonedma = movingaverageseasonal(detrended, per);
% deseasoned =  detrended - deseasonedc;
% figure(3)
% plot(deseasoned)
% title('Deseasoned prices time series')
% 
% % Autocorrelation of deseasoned time series using average seasonal component
% ac1 = autocorrelation(deseasoned, maxtau);
% autlim = zalpha/sqrt(n);
% figure(4)
% clf
% hold on
% for ii=1:maxtau
%     plot(ac1(ii+1,1)*[1 1],[0 ac1(ii+1,2)],'b','linewidth',1.5)
% end
% plot([0 maxtau+1],[0 0],'k','linewidth',1.5)
% plot([0 maxtau+1],autlim*[1 1],'--c','linewidth',1.5)
% plot([0 maxtau+1],-autlim*[1 1],'--c','linewidth',1.5)
% xlabel('\tau')
% ylabel('r(\tau)')
% title('Autocorrelation of deseasoned time series')
% % tittxt = sprintf('deseasoned time series by ave seasonal comp (%d), Ljung-Box test',per);
% % figure(5)
% 
% 
