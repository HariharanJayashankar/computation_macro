% Computes and analyzes aggregate dynamics
% Reproduces Figure 2, Table 3, Figure 3, Table 4, Table 8, Figure 5, and Figure 6
%
% Thomas Winberry, February 15th, 2018

% Declare structures/variables as global (to be referenced later)
global oo_ M_ options_ var_list_ nMeasure

% Set nMeasure
nMeasure		= 4;

%----------------------------------------------------------------
% Initialize with a run of Dynare
%----------------------------------------------------------------

clearvars -except nMeasure
close all
clc

cd('./Auxiliary Functions/Dynamics');

% Set parameters
setParameters_results;

% Compute grids
computeGrids;
computePolynomials;

% Save economic parameters
save economicParameters.mat ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ cchi
	
% Save approximation parameters
save approximationParameters.mat nProd nCapital nState nProdQuadrature nCapitalQuadrature nStateQuadrature ...
	nMeasureCoefficients nMeasure prodMin prodMax capitalMin capitalMax nShocks
	
% Save grids
save grids.mat vShocksGrid vShocksWeights mStateGrid mQuadratureGrid vQuadratureWeights mProdPrimeQuadrature

% Save polynomials
save polynomials.mat mStatePoly mStatePoly2 vStatePolySquared aProdPrimePoly mQuadraturePoly aProdPrimeQuadraturePoly

nMeasure

% Run Dynare
eval(['dynare dynamicModel.mod noclearall nograph -DnMeasure=' num2str(nMeasure) ' -Dperturbation_order=' num2str(1)])

cd('./Auxiliary Functions/Dynamics');

%----------------------------------------------------------------
% Compute forecasting regressions
%----------------------------------------------------------------

% Initialize
N                           			= 30;     
vSsigmaQRange               	= linspace(1e-4,.04,N);
tBurn                       			= 5000; 
T                           			= 20000;
vShocks                     		= randn(tBurn + T,2);
mAggregateTFPSeries        	= zeros(tBurn+T,N);
aMomentsSeries              	= zeros(nMeasureCoefficients,tBurn+T,N);
mMarginalUtilitySeries      	= zeros(tBurn+T,N);
mAggregateQSeries           = zeros(tBurn+T,N);



% Loop over parameters 
for iIndex = 1 : N

	% Set parameter value
	set_param_value('ssigmaQ',vSsigmaQRange(iIndex));
	
	% Solve the model
	info 												= stoch_simul(var_list_);
	
	% Simulate the model
	vYSeries 										= simult_(oo_.dr.ys,oo_.dr,vShocks,options_.order);
	
	% Extract series of interest
	mAggregateTFPSeries(:,iIndex)      = vYSeries(nState+nProd+2*nMeasureCoefficients+3,2:tBurn+T+1);
    mAggregateQSeries(:,iIndex)        	= vYSeries(nState+nProd+2*nMeasureCoefficients+4,2:tBurn+T+1);
	aMomentsSeries(:,:,iIndex)         	= vYSeries(nState+nProd+1:nState+nProd+nMeasureCoefficients,2:tBurn+T+1);
	mMarginalUtilitySeries(:,iIndex)   	= vYSeries(nState+nProd+2*nMeasureCoefficients+2,2:tBurn+T+1);
	
end


%%%
% Compute statistics of interest
%%%

% Preallocate matrices to store results
vR2Capital              		= zeros(N,1);
vR2MarginalUtility      	= zeros(N,1);
mKPrimeForecast         	= zeros(T,N);
mPForecast              		= zeros(T,N);
vDHMeanKPrime           	= zeros(N,1); 
vDHMaxKPrime            	= zeros(N,1);
vDHMeanP                		= zeros(N,1); 
vDHMaxP                 		= zeros(N,1);

% Loop over parameterizations
for iIndex = 1 : N
	
	% Compute aggregates
	vAggregateTFPSeries        		= mAggregateTFPSeries(:,iIndex);
	vAggregateQSeries          		= mAggregateQSeries(:,iIndex);
	mAggregateMomentsSeries    	= squeeze(aMomentsSeries(:,:,iIndex));
	vMarginalUtilitySeries     		= mMarginalUtilitySeries(:,iIndex);
	
	% Run simple regression for k_prime
	[b,bint,r,rint,stats]      			= regress(mAggregateMomentsSeries(2,tBurn+1:tBurn+T)',[ones(T,1),vAggregateTFPSeries(tBurn+1:tBurn+T),...
													vAggregateQSeries(tBurn+1:tBurn+T),mAggregateMomentsSeries(2,tBurn:tBurn+T-1)']);
	vR2Capital(iIndex,1)       		= stats(1);
	
	% Compute forecasted k_prime by forward iteration
	mKPrimeForecast(1,iIndex)  	= b(1) + b(2) * vAggregateTFPSeries(tBurn+1) + b(3) * vAggregateQSeries(tBurn+1) + ...
													b(4) * mAggregateMomentsSeries(2,tBurn);
	for t = 2 : T
		mKPrimeForecast(t,iIndex) = b(1) + b(2) * vAggregateTFPSeries(tBurn+t) + b(3) * vAggregateQSeries(tBurn+t) + ...
													b(4) * mKPrimeForecast(t-1,iIndex);
	end
	
	% Record DH Statistics
	vDHMeanKPrime(iIndex,1) 		= mean(abs(mKPrimeForecast(:,iIndex) - mAggregateMomentsSeries(2,tBurn+1:tBurn+T)'));
	vDHMaxKPrime(iIndex,1) 		= max(abs(mKPrimeForecast(:,iIndex) - mAggregateMomentsSeries(2,tBurn+1:tBurn+T)'));
	
	% Run simple regression for marginal utility
	[b,bint,r,rint,stats] 				= regress(log(vMarginalUtilitySeries(tBurn+1:tBurn+T)),[ones(T,1),vAggregateTFPSeries(tBurn+1:tBurn+T),...
													vAggregateQSeries(tBurn+1:tBurn+T),mAggregateMomentsSeries(2,tBurn:tBurn+T-1)']);
	vR2MarginalUtility(iIndex,1) 	= stats(1);
	
	% Compute forecasted p by forward iteration
	mPForecast(1,iIndex) 			= exp(b(1) + b(2) * vAggregateTFPSeries(tBurn+1) + b(3) * vAggregateQSeries(tBurn+t) + ...
													b(4) * mAggregateMomentsSeries(2,tBurn));
	for t = 2 : T
		mPForecast(t,iIndex) 			= exp(b(1) + b(2) * vAggregateTFPSeries(tBurn+t) + b(3) * vAggregateQSeries(tBurn+t) + ...
													b(4) * mKPrimeForecast(t-1,iIndex));
	end
	
	% Record DH Statistics
	vDHMeanP(iIndex,1) 				= mean(abs(mPForecast(:,iIndex) - vMarginalUtilitySeries(tBurn+1:tBurn+T))) ./ ...
													mean(vMarginalUtilitySeries(tBurn+1:tBurn+T));
	vDHMaxP(iIndex,1) 				= max(abs(mPForecast(:,iIndex) - vMarginalUtilitySeries(tBurn+1:tBurn+T))) ./ ...
													mean(vMarginalUtilitySeries(tBurn+1:tBurn+T));
		
end


%%%
% Plot
%%%

figure

subplot(1,3,1)
hold on
plot(vSsigmaQRange,vR2Capital,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vSsigmaQRange,vR2MarginalUtility,'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
h = legend('Aggregate capital','Marginal utility');
set(h,'interpreter','latex','location','southwest')
set(gcf,'color','w')
xlim([min(vSsigmaQRange) max(vSsigmaQRange)])
xlabel('Standard deviation $\sigma_{q}$','interpreter','latex')
ylabel('$R^{2}$','interpreter','latex')
title('$R^{2}$ of forecasting regressions','interpreter','latex','fontsize',14)
grid on
hold off

subplot(1,3,2)
hold on
plot(vSsigmaQRange,100*vDHMeanKPrime,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vSsigmaQRange,100*vDHMaxKPrime,'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
h = legend('Mean','Max');
set(h,'interpreter','latex','location','northwest')
set(gcf,'color','w')
xlim([min(vSsigmaQRange) max(vSsigmaQRange)])
xlabel('Standard deviation $\sigma_{q}$','interpreter','latex')
ylabel('Difference in mean log capital','interpreter','latex')
title('Maximum DH Error, Capital','interpreter','latex','fontsize',14)
grid on
hold off

subplot(1,3,3)
hold on
plot(vSsigmaQRange,100*vDHMeanP,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vSsigmaQRange,100*vDHMaxP,'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
set(gcf,'color','w')
xlim([min(vSsigmaQRange) max(vSsigmaQRange)])
xlabel('Standard deviation $\sigma_{q}$','interpreter','latex')
ylabel('Percent of mean marginal utility','interpreter','latex')
title('Maximum DH Error, Marginal Utility','interpreter','latex','fontsize',14)
grid on
hold off