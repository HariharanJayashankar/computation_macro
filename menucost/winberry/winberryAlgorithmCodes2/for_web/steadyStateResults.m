% Computes and analyzes steady state with no aggregate shocks
% Reproduces Figure 1 and Table 2 of paper
%
% Thomas Winberry, February 10th, 2018

% Maximum order of measure to consider
nOrderMax				= 6;


%----------------------------------------------------------------
% Compute steady state for different orders of approximation of distribution
%----------------------------------------------------------------


for iMeasureOrder	= 1 : nOrderMax

	clearvars -except iMeasureOrder nOrderMax
	close all
	clc
	
    % Print current order of measure
	iMeasureOrder

	cd('./Auxiliary Functions/Steady State');

	
	%----------------------------------------------------------------
	% Set parameters
	%----------------------------------------------------------------

	setParametersFigures;
	
	if iMeasureOrder == 1      % if first run, preallocate matrices to save results
	
		% Aggregates
		vOutput				= zeros(nOrderMax,1);
		vConsumption		= zeros(nOrderMax,1);
		vCapital			= zeros(nOrderMax,1);
		vWage				= zeros(nOrderMax,1);
		vMarginalUtility	= zeros(nOrderMax,1);
		save aggregates.mat vOutput vConsumption vCapital vWage vMarginalUtility

		% Distributions
		aFineDistribution	= zeros(nProdFine,nCapitalFine,nOrderMax);
		save distributions.mat aFineDistribution
		
	end

	
	%----------------------------------------------------------------
	% Solve for steady state wage
	%----------------------------------------------------------------

	coreSteadyState;

	
	%----------------------------------------------------------------
	% Compute steady state objects over histogram grid for plots
	%----------------------------------------------------------------

	[~,vCoefficients,vCapitalCoefficients,vCapitalAdjust,...
		vCapitalConstrained,vCutoff,vHistogram] = computeLMCResidualHistogram(wage);
		
		
	%----------------------------------------------------------------
	% Compute aggregates variables
	%----------------------------------------------------------------
	[resid,vMoments,vParameters,aggregateConsumption,marginalUtility,aggregateOutput,aggregateCapital,aggregateInvestment] = ...
		computeLMCResidualPolynomials(wage,vParameters,vMomentsHistogram,mGridMoments);
	aggregateHours = resid + nSS;
		
		
	%----------------------------------------------------------------
	% Compute density from parametric family over fine grid
	%----------------------------------------------------------------

	% Compute fine grid of centered moments
	mGridMomentsFine 			= zeros(nStateFine,2);
	mGridMomentsFine(:,1) 	= mFineGrid(:,1) - vMoments(1,1);
	mGridMomentsFine(:,2) 	= log(mFineGrid(:,2)) - vMoments(2,1);

	counter = 3;
	for i = 2:nMeasure
		for j = 0:i
			mGridMomentsFine 	= [mGridMomentsFine ((mFineGrid(:,1) - vMoments(1,1)) .^ (i - j)) .* ...
											((log(mFineGrid(:,2)) - vMoments(2,1)) .^ j) - vMoments(counter,1)];
			counter 					= counter + 1;
		end
	end

	% Compute distribution over fine grid
	vFineDistribution 			= vParameters(1,1) * exp(mGridMomentsFine * vParameters(2:nMeasureCoefficients+1));
	
	
	%----------------------------------------------------------------
	% Save output from this run
	%----------------------------------------------------------------
	
	% Aggregates
	load aggregates.mat
	vOutput(iMeasureOrder,1)			= aggregateOutput;
	vConsumption(iMeasureOrder,1)	= aggregateConsumption;
	vCapital(iMeasureOrder,1)			= aggregateCapital;
	vWage(iMeasureOrder,1)				= wage;
	vMarginalUtility(iMeasureOrder,1)	= marginalUtility;
	save aggregates.mat vOutput vConsumption vCapital vWage vMarginalUtility
	
	% Distributions
	load distributions.mat
	aFineDistribution(:,:,iMeasureOrder) = reshape(vFineDistribution,nProdFine,nCapitalFine);
	save distributions.mat aFineDistribution
	
	cd('../../')
	
end

cd('./Auxiliary Functions/Steady State');


%----------------------------------------------------------------
% Compute steady state using histogram (for comparison of aggregates)
%----------------------------------------------------------------

% Solve for market clearing wage using histogram
f 			= @(wage) computeLMCResidualHistogram(wage);
options 	= optimoptions('fsolve','Display','iter-detailed');
[wageHistogram,err,exitflag] = fsolve(f,wage,options);

% Compute steady state objects under histogram approximation
[~,vCoefficients,vCapitalCoefficients,vCapitalAdjust,...
	vCapitalConstrained,vCutoff,vHistogram] = computeLMCResidualHistogram(wageHistogram);
[vCapitalAdjust,vCapitalConstrained,vCutoff] = ...
	computePolicies(vCoefficients,vCapitalCoefficients,wageHistogram,mFineGrid,mFinePoly,aProdPrimeFinePoly);
mHistogram = reshape(vHistogram,nProdFine,nCapitalFine);

% Aggregate consumption 
vLaborDemand 	= ((exp(mFineGrid(:,1)) .* (mFineGrid(:,2) .^ ttheta) * nnu) / wageHistogram) .^ (1 / (1 - nnu));
vIntegrand 		= exp(mFineGrid(:,1)) .* (mFineGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu) + ...
							(vCutoff ./ ppsiCapital) .* (-(vCapitalAdjust - (1 - ddelta) * mFineGrid(:,2))) + (1 - (vCutoff ./ ppsiCapital)) .* ...
							(-(vCapitalConstrained - (1 - ddelta) * mFineGrid(:,2)));
aggregateConsumptionHistogram = sum(vIntegrand(:) .* vHistogram(:));

% Marginal utility
marginalUtilityHistogram 	= aggregateConsumption ^ (-ssigma);

% Output
vIntegrand 						= exp(mFineGrid(:,1)) .* (mFineGrid(:,2) .^ ttheta) .* (vLaborDemand .^ nnu);
aggregateOutputHistogram	= sum(vIntegrand(:) .* vHistogram(:));

% Capital
aggregateCapitalHistogram	= sum(mFineGrid(:,2) .* vHistogram(:));

% Investment
vIntegrand 						= (vCutoff ./ ppsiCapital) .* (vCapitalAdjust - (1 - ddelta) * mFineGrid(:,2)) + (1 - (vCutoff ./ ppsiCapital)) .* ...
											(vCapitalConstrained - (1 - ddelta) * mFineGrid(:,2));
aggregateInvestmentHistogram = sum(vIntegrand(:) .* vHistogram(:));

%----------------------------------------------------------------
% Plot and analyze steady state objects
%----------------------------------------------------------------

%%%
% Print out aggregates (Table 2)
%%%

% Parametric family
load aggregates.mat
vOutput				
vConsumption	
vCapital		
vWage				
vMarginalUtility	

% Histogram
aggregateOutputHistogram
aggregateConsumptionHistogram
aggregateCapitalHistogram
wageHistogram
marginalUtilityHistogram


%%%
% Comparison of Marginal Distributions
%%%

load distributions.mat

figure

subplot(2,2,1)
hold on
plot(vProdGridFine,sum(squeeze(aFineDistribution(:,:,1)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,1)))),...
	'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
h1 = scatter(vProdGridFine(1:nProdFine/5:nProdFine),sum(squeeze(aFineDistribution(1:nProdFine/5:nProdFine,:,1)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,1)))),...
	'x','MarkerFaceColor',[178/255,34/255,34/255],'MarkerEdgeColor',[178/255,34/255,34/255]);
plot(vProdGridFine,sum(squeeze(aFineDistribution(:,:,2)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,2)))),...
	'linewidth',1.5,'linestyle','--','color',[20*8/255,62/255,.2*118/255])
h2 = scatter(vProdGridFine(1:nProdFine/5:nProdFine),sum(squeeze(aFineDistribution(1:nProdFine/5:nProdFine,:,2)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,2)))),...
	'+','MarkerFaceColor',[20*8/255,62/255,.2*118/255],'MarkerEdgeColor',[20*8/255,62/255,.2*118/255]);
plot(vProdGridFine,sum(squeeze(aFineDistribution(:,:,4)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,4)))),...
	'linewidth',1.5,'linestyle','--','color',[10*8/255,62/255,.6*118/255])
h3 = scatter(vProdGridFine(1:nProdFine/5:nProdFine),sum(squeeze(aFineDistribution(1:nProdFine/5:nProdFine,:,4)),2) ./ sum(sum(squeeze(aFineDistribution(:,:,4)))),...
	'o','MarkerFaceColor',[10*8/255,62/255,.6*118/255],'MarkerEdgeColor',[10*8/255,62/255,.6*118/255]);
h4 = plot(vProdGridFine,sum(mHistogram,2),'linewidth',1.5,'color',[8/255,62/255,118/255],'linestyle','-');
grid on
xlabel('Productivity, $\varepsilon$','interpreter','latex')
ylabel('Density of firms, $g(\varepsilon)$','interpreter','latex')
xlim([.9*prodMin .9*prodMax])
title('Marginal distribution of productivity','interpreter','latex','fontsize',12)
set(gcf,'color','w')
hold off

subplot(2,2,2)
hold on
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(:,:,1)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,1)))),...
	'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
h1 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),sum(squeeze(aFineDistribution(:,1:nCapitalFine/5:nCapitalFine,1)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,1)))),...
	'x','MarkerFaceColor',[178/255,34/255,34/255],'MarkerEdgeColor',[178/255,34/255,34/255]);
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(:,:,2)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,2)))),...
	'linewidth',1.5,'linestyle','--','color',[20*8/255,62/255,.2*118/255])
h2 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),sum(squeeze(aFineDistribution(:,1:nCapitalFine/5:nCapitalFine,2)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,2)))),...
	'+','MarkerFaceColor',[20*8/255,62/255,.2*118/255],'MarkerEdgeColor',[20*8/255,62/255,.2*118/255]);
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(:,:,4)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,4)))),...
	'linewidth',1.5,'linestyle','--','color',[10*8/255,62/255,.6*118/255])
h3 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),sum(squeeze(aFineDistribution(:,1:nCapitalFine/5:nCapitalFine,4)),1) ./ sum(sum(squeeze(aFineDistribution(:,:,4)))),...
	'o','MarkerFaceColor',[10*8/255,62/255,.6*118/255],'MarkerEdgeColor',[10*8/255,62/255,.6*118/255]);
h4 = plot(vCapitalGridFine,sum(mHistogram,1),'linewidth',1.5,'color',[8/255,62/255,118/255],'linestyle','-');
h = legend([h1,h2,h3,h4],'$n_{g} = 1$','$n_{g} = 2$','$n_{g} = 4$','Histogram');
set(h,'interpreter','latex','fontsize',12,'location','northeast')
set(gcf,'color','w')
grid on
xlabel('Capital, $k$','interpreter','latex')
ylabel('Density of firms, $g(\varepsilon)$','interpreter','latex')
title('Marginal distribution of capital','interpreter','latex','fontsize',12)
xlim([1.1*capitalMin .9*capitalMax])
hold off

subplot(2,2,3)
hold on
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,2)),1) ./ sum(sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,2)))),...
	'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
h1 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),squeeze(aFineDistribution(ceil(nProdFine/4),1:nCapitalFine/5:nCapitalFine,2)) ./ sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,2))),...
	'x','MarkerFaceColor',[178/255,34/255,34/255],'MarkerEdgeColor',[178/255,34/255,34/255]);
h2 = plot(vCapitalGridFine,mHistogram(ceil(nProdFine/4),:) / sum(mHistogram(ceil(nProdFine/4),:)),'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','-');
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,2)),1) ./ sum(sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,2)))),...
	'linewidth',1.5,'linestyle','--','color',[8/255,62/255,118/255])
h3 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),squeeze(aFineDistribution(ceil(3*nProdFine/4),1:nCapitalFine/5:nCapitalFine,2)) ./ sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,2))),...
	'x','MarkerFaceColor',[8/255,62/255,118/255],'MarkerEdgeColor',[8/255,62/255,118/255]);
h4 = plot(vCapitalGridFine,mHistogram(ceil(3*nProdFine/4),:) / sum(mHistogram(ceil(3*nProdFine/4),:)),'linewidth',1.5,'color',[8/255,62/255,118/255],'linestyle','-');
h = legend([h2,h4],'Low productivity','High productivity');
set(h,'interpreter','latex','fontsize',12,'location','northeast')
set(gcf,'color','w')
grid on
xlabel('Capital, $k$','interpreter','latex')
ylabel('Density of firms, $g(k|\varepsilon)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Conditional distribution of capital, $n_{g} = 2$','interpreter','latex','fontsize',12)
hold off

subplot(2,2,4)
hold on
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,4)),1) ./ sum(sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,4)))),...
	'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
h1 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),squeeze(aFineDistribution(ceil(nProdFine/4),1:nCapitalFine/5:nCapitalFine,4)) ./ sum(squeeze(aFineDistribution(ceil(nProdFine/4),:,4))),...
	'x','MarkerFaceColor',[178/255,34/255,34/255],'MarkerEdgeColor',[178/255,34/255,34/255]);
h2 = plot(vCapitalGridFine,mHistogram(ceil(nProdFine/4),:) / sum(mHistogram(ceil(nProdFine/4),:)),'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','-');
plot(vCapitalGridFine,sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,4)),1) ./ sum(sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,4)))),...
	'linewidth',1.5,'linestyle','--','color',[8/255,62/255,118/255])
h3 = scatter(vCapitalGridFine(1:nCapitalFine/5:nCapitalFine),squeeze(aFineDistribution(ceil(3*nProdFine/4),1:nCapitalFine/5:nCapitalFine,4)) ./ sum(squeeze(aFineDistribution(ceil(3*nProdFine/4),:,4))),...
	'x','MarkerFaceColor',[8/255,62/255,118/255],'MarkerEdgeColor',[8/255,62/255,118/255]);
h4 = plot(vCapitalGridFine,mHistogram(ceil(3*nProdFine/4),:) / sum(mHistogram(ceil(3*nProdFine/4),:)),'linewidth',1.5,'color',[8/255,62/255,118/255],'linestyle','-');
set(gcf,'color','w')
grid on
xlabel('Capital, $k$','interpreter','latex')
ylabel('Density of firms, $g(k|\varepsilon)$','interpreter','latex')
title('Conditional distribution of capital, $n_{g} = 4$','interpreter','latex','fontsize',12)
xlim([1.1*capitalMin .9*capitalMax])
hold off



cd('../../')
