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
figure(8);
plot(nrmseV,"-o");
title("NRMSE");
xlabel("steps forward");
ylabel("value");


%% LINEAR 4
% %predictions=[];
% indexes=[];
% nrmse=[];
% step=30;
% for n1 = 20:step:n %6 different points in total
%    indexes=[indexes n1];
%    
%    [nrmseV,~,~,~] = predictARMAnrmse(deseasoned,5,0,1,n-n1,[]);
%    nrmse=[nrmse nrmseV];
%    
%    %[predictions] =[predictions predictARMAmultistep(deseasoned,n1,5,0,1,[])];
%    %[preV] = predictARMAmultistep(xV,n1,p,q,Tmax,'example');
% end
% 
% figure()
% 
% % plot(indexes,predictions,"-o");
% % hold on;
% plot(indexes,nrmse,"-o");
% title("NRMSE for one step forward prediction with AR(5) model, in prices time series");
% xlabel("Size of learning test");
% legend("NRMSE");

index=1;
nrmseAR5 = zeros(5,1);
for i=70:5:90
    [nrmseAR5(index),~,~,~] = predictARMAnrmse(deseasoned, 5, 0, 1, round(0.01*(100-i)*n), '');
    index=index+1;
end
figure(9)
plot((0.7:0.05:0.9).',nrmseAR5,'-o')
title('NRMSE for one step forward price prediction with AR(5) model')
ylabel('NRMSE value')
xlabel('Percentage of time series used as training set ')

%% NO-LINEAR 
%Getting the residuals from our best linear model fit( in linear 3).
residuals =deseasoned - predict(armamodel,deseasoned)
figure(10);
plot(residuals);
title('Residuals from the linear model fit( in linear 3). ');
xlabel('t');
ylabel('residual');

iidtimeSeries=residuals

%Resampling the residuals time series
for i=1:20 
    RandIndexes=randperm(n-1);
    iidtimeSeries=[iidtimeSeries  residuals(RandIndexes)];
end

%Linear autocorrelation
%maxtau1 for the autocorrelation and mutual information
maxtau1=3;
nOfnewSamples=20;
LinearAutoCorrelations=zeros(maxtau1,nOfnewSamples+1);

%computing autocorrelation for all the new samples
figure()
for i=1:(nOfnewSamples+1) 
    autoCorr=autocorrelation(iidtimeSeries(:,i), maxtau1)
    LinearAutoCorrelations(:,i)=autoCorr(2:(maxtau1+1),2);
end
clf;

%plotting the autocorrelation function
figure(10)
for i=1:maxtau1
    plot([0:nOfnewSamples] , LinearAutoCorrelations(i,:),'-o');
    hold on
end
plot([0 nOfnewSamples+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 nOfnewSamples+1],-autlim*[1 1],'--c','linewidth',1.5)
legend("\tau=1","\tau=2","\tau=3")
title('Linear autocorrelation for the 21 samples.')
ylabel("autocorrelation")
xlabel("sample Number (n=0 belongs to the residuals)");


%Mutual information for the 21 samples
SamplesMutualInfo=zeros(maxtau1,nOfnewSamples+1);
figure()
%computing mutual information for all the new samples
for i=1:(nOfnewSamples+1)
    mut=mutualinformation(iidtimeSeries(:,i), maxtau1)
    SamplesMutualInfo(:,i)=mut(2:(maxtau1+1),2);
end
clf;

%plotting the mutual information function
figure(11)
for i=1:maxtau1
    plot([0:nOfnewSamples] , SamplesMutualInfo(i,:),'-o');
    hold on
end
plot([0 nOfnewSamples+1],autlim*[1 1],'--c','linewidth',1.5)
legend("\tau=1","\tau=2","\tau=3")
title('Mutual information for the 21 samples.')
ylabel("mutual information")
xlabel("sample Number (0 belongs to the original residuals)");

%Histogram for lat(tau) =1 for autocorrelation and mutual information
% bins=10;
figure(12)
histogram(LinearAutoCorrelations(1,2:(nOfnewSamples+1)));
hold on;
line([LinearAutoCorrelations(1,1) LinearAutoCorrelations(1,1)], [0 9],'Color','red','linewidth',1.5);
title("Autocorrealtion for the 20 new samples and for original residuals of price(red)");

figure(13)
histogram(SamplesMutualInfo(1,2:(nOfnewSamples+1)));
hold on;
line([SamplesMutualInfo(1,1) SamplesMutualInfo(1,1)], [0 9],'Color','red','linewidth',1.5);
title("Mutual information for the 20 new samples and for original residuals of price(red)");