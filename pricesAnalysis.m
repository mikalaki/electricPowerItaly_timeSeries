%% Ergasia Xronoseirwn : Zisou Charilaos AEM 9213 ,Karatzas Michalis AEM 9137
%% Linear Analysis
close all; clear; clc;

% import data
load('pricesData.mat');

% %getting the time series (START)
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
%  %getting the time series (END)


%constants
per = 7;
maxtau = 100;
alpha = 0.05;
Tmax=20;
zalpha = norminv(1-alpha/2);
n = length(prices)-1;
bins=10

% plot unprocessed time series
figure()
plot(prices)
xlabel('t')
ylabel('y(t)')
title('Unprocessed prices time series')

%% LINEAR 1 
% 1a. Detrend time series by first order differences
detrended=zeros(n,1);
for i=1:n
    detrended(i)=prices(i+1)-prices(i);
end
figure()
plot(detrended)
title('Detrended prices time series by first order differences')

%1b. Deseason time series
deseasonedc = seasonalcomponents(detrended, per);
deseasoned =  detrended - deseasonedc;
figure()
plot(deseasoned)
title('Deseasoned prices time series')


%%  LINEAR 2
%2. Autocorrelation of deseasoned time series 
ac1 = autocorrelation(deseasoned(~isnan(deseasoned)), maxtau);
autlim = zalpha/sqrt(n);
figure()
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
title('Autocorrelation of prices time series')
hold off

% Partial autocorrelation
pac1 = parautocor(deseasoned, maxtau);
figure()
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
title('Partial autocorrelation of prices time series')
hold off

% Ljung-Box Portmanteau Test
figure()
tittxt = ('Prices time series');
[h2V,p2V,Q2V] = portmanteauLB(ac1(2:maxtau+1,2),n,alpha,tittxt);



%% LINEAR 3
% AIC
%aicMatrix = akaike(deseasoned, Tmax);


%ARMA model parameters.
p=3
q=3

%fitting the model
[nrmseARMA,phiV,thetaV,SDz,aicS,fpeS,armamodel]=fitARMA(deseasoned, p, q, Tmax);

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
plot(nrmseARMA,"-o");
title('Prices: NRMSE of ARMA(3,3) for T=1,2...10')
xlabel("T");
ylabel("NRMSE");


%% LINEAR 4
index=1;
nrmseAR5 = zeros(5,1);
for i=70:5:90
    [nrmseAR5(index),~,~,~] = predictARMAnrmse(deseasoned, 5, 0, 1, round(0.01*(100-i)*n), '');
    index=index+1;
end
figure()
plot((0.7:0.05:0.9).',nrmseAR5,'-o')
ylabel("NRMSE");
xlabel("percent of training data");
title('Prices: NRMSE of AR(5) using 70%, 75%, 80%, 85%, 90% of training data')

%% Non-Linear Analysis
%Getting the residuals from our best linear model fit (in linear 3)
residuals = deseasoned - predict(armamodel, deseasoned);
resSamples=residuals;

% Resampling the residuals time series
nOfnewSamples=20;
for i=1:nOfnewSamples 
    RandIndexes=randperm(n);
    resSamples=[resSamples  residuals(RandIndexes)];
end

%Linear autocorrelation
%maxtau1 for the autocorrelation and mutual information plot
maxtau1=4;
LinearAutoCorrelations=zeros(maxtau1,nOfnewSamples+1);

%computing autocorrelation for all the new samples and original residues
for i=1:(nOfnewSamples+1) 
    autoCorr=autocorrelation(resSamples(:,i), maxtau1)
    LinearAutoCorrelations(:,i)=autoCorr(2:(maxtau1+1),2);
end


%plotting the autocorrelation function
figure()
for i=1:maxtau1
    plot([0:nOfnewSamples] , LinearAutoCorrelations(i,:),'-o');
    hold on
end
plot([0 nOfnewSamples+1],autlim*[1 1],'--c','linewidth',1.5)
plot([0 nOfnewSamples+1],-autlim*[1 1],'--c','linewidth',1.5)
legend("\tau=1","\tau=2","\tau=3","\tau=4")
title('Linear autocorrelation for the 21 samples.(price)')
ylabel("autocorrelation")
xlabel("Sample number (0 belongs to the original residuals, 1-20 to permutations)");


%Mutual information for the 21 samples(original +permutations)
SamplesMutualInfo=zeros(maxtau1,nOfnewSamples+1);

figure()
%computing mutual information for all the new samples
for i=1:(nOfnewSamples+1)
    mut=mutualinformation(resSamples(:,i), maxtau1)
    SamplesMutualInfo(:,i)=mut(2:(maxtau1+1),2);
end


%plotting the mutual information function
figure()
for i=1:maxtau1
    plot([0:nOfnewSamples] , SamplesMutualInfo(i,:),'-o');
    hold on
end

legend("\tau=1","\tau=2","\tau=3","\tau=4")
title('Mutual information for the 21 samples. (price)')
ylabel("mutual information")
xlabel("Sample number (0 belongs to the original residuals, 1-20 to permutations)");

%Histogram for lat(tau) =1 for autocorrelation 
figure()
histogram(LinearAutoCorrelations(1,2:(nOfnewSamples+1)),bins);
hold on;
line([LinearAutoCorrelations(1,1) LinearAutoCorrelations(1,1)], [0 9],'Color','red','linewidth',1.5);
title("AC histogram of new samples and original residues (red) (price)");
ylabel("frequency");
xlabel("Autocorrealtion value");

%Histogram for lat(tau) =1 for Mutual Information
figure()
histogram(SamplesMutualInfo(1,2:(nOfnewSamples+1)),bins);
hold on;
line([SamplesMutualInfo(1,1) SamplesMutualInfo(1,1)], [0 9],'Color','red','linewidth',1.5);
title("MI histogram of new samples and original residues(red) (price)");
ylabel("frequency");
xlabel("Mutual Information value");

%lag
tau=3;

%FNN
figure();
fnn = falsenearest(residuals,tau,10,10,0,'Price Residuals');

%getting correlation dimension
SamplesCorrDimension=zeros(1,nOfnewSamples+1);

%embedding dimensions
m=4;

%correlation dimension for the residuals
figure();
[~,~,~,~,~] = correlationdimension(residuals,tau,10,"Prices residuals");
 
 %%with log(r1), logf(r2) and fac fixed on the series.
 %figure();
 %[~,~,~,~,~] = correlationdimension(residuals,1,10,"Prices residuals",0.2,1.1,1.7);

%correlation dimension for the 20 samples ,we get embeding distance m=4



for i=1:(nOfnewSamples+1)
    [~,~,~,~,nuM] = correlationdimension_no_plot(resSamples(:,i),tau,10,"Prices residuals");
    SamplesCorrDimension(:,i) = nuM(m,4);
end

%Histogram for correlation dimension, with embedding dimension m=4 and lag τ=4
figure()
histogram(SamplesCorrDimension(1,2:(nOfnewSamples+1)),bins);
hold on;
line([SamplesCorrDimension(1,1) SamplesCorrDimension(1,1)], [0 9],'Color','red','linewidth',1.5);
title("CD histogram of new samples and original residues(red) (price)");
ylabel("frequency");
xlabel("Correlation Dimension value");



