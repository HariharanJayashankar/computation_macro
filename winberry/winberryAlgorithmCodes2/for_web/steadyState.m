% Computes and analyzes steady state with no aggregate shocks.
% 
%
% Thomas Winberry, February 14th, 2018

% Clean up workspace
clear all
close all
clc

% Change directory 
cd('./Auxiliary Functions/Steady State');


%----------------------------------------------------------------
% Set parameters
%----------------------------------------------------------------

setParameters;


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
% Compute marginal value function 
%----------------------------------------------------------------

% Compute adjust capital accumulation policy
[vCapitalAdjust,vCapitalConstrained,vCutoff] = computePolicies(vCoefficients,vCapitalCoefficients,wage,mStateGrid,mStatePoly,aProdPrimePoly);

% Compute coefficients
mCoefficientsDeriv = reshape(vCoefficients,nProd,nCapital);
mCoefficientsDeriv = (2 / (capitalMax - capitalMin)) * repmat(linspace(1,nCapital-1,nCapital-1),[nProd 1]) .* ...
	mCoefficientsDeriv(:,2:nCapital);
vCoefficientsDeriv = reshape(mCoefficientsDeriv,nState-nProd,1);


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
mGridMomentsFine		 = zeros(nStateFine,2);
mGridMomentsFine(:,1) = mFineGrid(:,1) - vMoments(1,1);
mGridMomentsFine(:,2) = log(mFineGrid(:,2)) - vMoments(2,1);

counter = 3;
for i = 2:nMeasure
	for j = 0:i
		mGridMomentsFine 	= [mGridMomentsFine ((mFineGrid(:,1) - vMoments(1,1)) .^ (i - j)) .* ...
			((log(mFineGrid(:,2)) - vMoments(2,1)) .^ j) - vMoments(counter,1)];
		counter 					= counter + 1;
	end
end

% Compute distribution over fine grid
vFineDistribution 	= vParameters(1,1) * exp(mGridMomentsFine * vParameters(2:nMeasureCoefficients+1));
mFineDistribution 	= reshape(vFineDistribution,nProdFine,nCapitalFine);


%----------------------------------------------------------------
% Plot and analyze steady state objects
%----------------------------------------------------------------

% Compute other objects over fine grid for plotting
[vCapitalAdjust,vCapitalConstrained,vCutoff] 	= computePolicies(vCoefficients,vCapitalCoefficients,wage,mFineGrid,mFinePoly,aProdPrimeFinePoly);
mValueFunction 											= reshape(mFinePoly * vCoefficients,nProdFine,nCapitalFine);
mMarginalValueFunction 								= reshape(mFinePoly2 * vCoefficientsDeriv,nProdFine,nCapitalFine);
mCapitalAdjust 											= reshape(vCapitalAdjust,nProdFine,nCapitalFine);
mCutoff 													= reshape(vCutoff,nProdFine,nCapitalFine);
mHistogram 												= reshape(vHistogram,nProdFine,nCapitalFine);


% Print steady state aggregates
fprintf('----- Steady state aggregates ------- \n')
fprintf('Aggregate Output: %5.3f\n',aggregateOutput)
fprintf('Aggregate Consumption: %5.3f\n',aggregateConsumption)
fprintf('Aggregate Investment: %5.3f\n',aggregateInvestment)
fprintf('Aggregate Capital: %5.3f\n',aggregateCapital)
fprintf('Wage: %5.3f\n',wage)
fprintf('Marginal utility: %5.3f\n\n',marginalUtility)

% Plot properties of individual decisions
figure

subplot(2,2,1)
hold on
plot(vCapitalGridFine,mValueFunction(ceil(nProdFine/4),:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mValueFunction(ceil(3*nProdFine/4),:),'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Capital, $k$','interpreter','latex')
ylabel('Value function, $v(\varepsilon,k)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
h	 = legend('Low productivity','High productivity','location','southeast');
set(h,'interpreter','latex','location','northwest','fontsize',12)
set(gcf,'color','w')
title('Value function','interpreter','latex','fontsize',14)
grid on
hold off

subplot(2,2,2)
hold on
plot(vCapitalGridFine,mMarginalValueFunction(ceil(nProdFine/4),:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mMarginalValueFunction(ceil(3*nProdFine/4),:),'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Capital, $k$','interpreter','latex')
ylabel('Marginal value function, $\partial v(\varepsilon,k) / \partial k $','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
set(gcf,'color','w')
title('Marginal value of capital','interpreter','latex','fontsize',14)
grid on
hold off

subplot(2,2,3)
hold on
plot(vCapitalGridFine,mCapitalAdjust(ceil(nProdFine/4),:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mCapitalAdjust(ceil(3*nProdFine/4),:),'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vCapitalGridFine,vCapitalGridFine,'k--','linewidth',1)
xlabel('Capital, $k$','interpreter','latex')
ylabel('Next period capital, $k(\varepsilon,k)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Investment policy, conditional on adjusting','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

subplot(2,2,4)
hold on
plot(vCapitalGridFine,mCutoff(ceil(nProdFine/4),:) ./ ppsiCapital,'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mCutoff(ceil(3*nProdFine/4),:) ./ ppsiCapital,'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Capital, $k$','interpreter','latex')
ylabel('Fixed cost threshold, $\xi(\varepsilon,k)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Probability of adjusting capital','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

% Plot distributions
figure

subplot(2,2,1)
hold on
plot(vProdGridFine,sum(mHistogram,2)','linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vProdGridFine,sum(mFineDistribution,2)' ./ sum(mFineDistribution(:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Productivity, $\varepsilon$','interpreter','latex')
ylabel('Mass of firms, $\mu(\varepsilon)$','interpreter','latex')
xlim([.9*prodMin .9*prodMax])
title('Invariant distribution of firms','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off

subplot(2,2,2)
hold on
plot(vCapitalGridFine,sum(mHistogram,1),'linewidth',1.5,'color',[8/255,62/255,118/255])
plot(vCapitalGridFine,sum(mFineDistribution,1) ./ sum(mFineDistribution(:)),...
	'linewidth',1.5,'color',[178/255,34/255,34/255],'linestyle','--')
xlabel('Capital, $k$','interpreter','latex')
ylabel('Mass of firms, $\mu(k)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Invariant distribution of firms','interpreter','latex','fontsize',14)
h	= legend('Histogram','Parametric Family');
set(h,'interpreter','latex','fontsize',12,'location','northeast')
set(gcf,'color','w')
grid on
hold off

subplot(2,2,3)
hold on
plot(vCapitalGridFine,mHistogram(ceil(nProdFine/4),:),...
	'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mHistogram(ceil(3*nProdFine/4),:),...
	'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Capital, $k$','interpreter','latex')
ylabel('Conditional mass of firms, $\mu(k|\varepsilon)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Conditional distribution, histogram','interpreter','latex','fontsize',14)
h	= legend('Low productivity','High productivity','location','northeast');
set(h,'interpreter','latex','location','northeast','fontsize',12)
set(gcf,'color','w')
grid on
hold off

subplot(2,2,4)
hold on
plot(vCapitalGridFine,mFineDistribution(ceil(nProdFine/4),:),'linewidth',1.5,'color',[178/255,34/255,34/255])
plot(vCapitalGridFine,mFineDistribution(ceil(3*nProdFine/4),:),'linewidth',1.5,'color',[8/255,62/255,118/255])
xlabel('Capital, $k$','interpreter','latex')
ylabel('Conditional mass of firms, $\mu(k|\varepsilon)$','interpreter','latex')
xlim([1.1*capitalMin .9*capitalMax])
title('Conditional distribution, parametric family','interpreter','latex','fontsize',14)
set(gcf,'color','w')
grid on
hold off



cd('../../')
