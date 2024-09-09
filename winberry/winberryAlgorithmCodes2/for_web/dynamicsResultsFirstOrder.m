% Computes and analyzes aggregate dynamics
% Reproduces Figure 2, Table 3, Figure 3, Table 4, Table 8, Figure 5, and Figure 6
%
% Thomas Winberry, February 15th, 2018

% Declare structures/variables as global (to be referenced later)
global oo_ M_ options_ var_list_ nMeasure


%----------------------------------------------------------------
% Loop over different degrees of approximation for distribution
%----------------------------------------------------------------

for nMeasure		= 2 : 4

	clearvars -except nMeasure
	close all
	
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
	eval(['dynare dynamicModel.mod noclearall -DnMeasure=' num2str(nMeasure) ' -Dperturbation_order=' num2str(1)])

	% Save aggregate IRFs output (have to loop over nMeasure manually for now)
	eval(sprintf('vOutput_TFP_Order_%d = oo_.irfs.logAggregateOutput_aggregateTFPShock;',nMeasure));
	eval(sprintf('vInvestment_TFP_Order_%d = oo_.irfs.logAggregateInvestment_aggregateTFPShock;',nMeasure));
	eval(sprintf('vConsumption_TFP_Order_%d = oo_.irfs.logAggregateConsumption_aggregateTFPShock;',nMeasure));
	eval(sprintf('vHours_TFP_Order_%d = oo_.irfs.logAggregateHours_aggregateTFPShock;',nMeasure));
	eval(sprintf('vWage_TFP_Order_%d = oo_.irfs.logWage_aggregateTFPShock;',nMeasure));
	eval(sprintf('vInterestRate_TFP_Order_%d = oo_.irfs.realInterestRate_aggregateTFPShock;',nMeasure));

	eval(sprintf('vOutput_Q_Order_%d = oo_.irfs.logAggregateOutput_aggregateQShock;',nMeasure));
	eval(sprintf('vInvestment_Q_Order_%d = oo_.irfs.logAggregateInvestment_aggregateQShock;',nMeasure));
	eval(sprintf('vConsumption_Q_Order_%d = oo_.irfs.logAggregateConsumption_aggregateQShock;',nMeasure));
	eval(sprintf('vHours_Q_Order_%d = oo_.irfs.logAggregateHours_aggregateQShock;',nMeasure));
	eval(sprintf('vWage_Q_Order_%d = oo_.irfs.logWage_aggregateQShock;',nMeasure));
	eval(sprintf('vInterestRate_Q_Order_%d = oo_.irfs.realInterestRate_aggregateQShock;',nMeasure));

	eval(sprintf('save first_order_irfs_%d.mat vOutput_TFP_Order_%d vInvestment_TFP_Order_%d vConsumption_TFP_Order_%d vHours_TFP_Order_%d vWage_TFP_Order_%d vInterestRate_TFP_Order_%d vOutput_Q_Order_%d vInvestment_Q_Order_%d vConsumption_Q_Order_%d vHours_Q_Order_%d vWage_Q_Order_%d vInterestRate_Q_Order_%d',nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure));

	% Save distributional IRFs
	eval(sprintf('vMoment_2_TFP_Order_%d = oo_.irfs.moment_2_aggregateTFPShock;',nMeasure));
	eval(sprintf('vMoment_4_TFP_Order_%d = oo_.irfs.moment_4_aggregateTFPShock;',nMeasure));
	eval(sprintf('vMoment_5_TFP_Order_%d = oo_.irfs.moment_5_aggregateTFPShock;',nMeasure));
	eval(sprintf('vLog_marginal_utility_TFP_Order_%d = oo_.irfs.logMarginalUtility_aggregateTFPShock;',nMeasure));

	eval(sprintf('vMoment_2_Q_Order_%d = oo_.irfs.moment_2_aggregateQShock;',nMeasure));
	eval(sprintf('vMoment_4_Q_Order_%d = oo_.irfs.moment_4_aggregateQShock;',nMeasure));
	eval(sprintf('vMoment_5_Q_Order_%d = oo_.irfs.moment_5_aggregateQShock;',nMeasure));
	eval(sprintf('vLog_marginal_utility_Q_Order_%d = oo_.irfs.logMarginalUtility_aggregateQShock;',nMeasure));

	eval(sprintf('save distributional_irfs_%d.mat vMoment_2_TFP_Order_%d vMoment_4_TFP_Order_%d vMoment_5_TFP_Order_%d vLog_marginal_utility_TFP_Order_%d vMoment_2_Q_Order_%d vMoment_4_Q_Order_%d vMoment_5_Q_Order_%d vLog_marginal_utility_Q_Order_%d',nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure,nMeasure));

end

return


%----------------------------------------------------------------
% Plot impulse responses for different orders of approximation
%----------------------------------------------------------------

% Load impulse responses
load first_order_irfs_2.mat
load first_order_irfs_3.mat
load first_order_irfs_4.mat
load distributional_irfs_2.mat
load distributional_irfs_3.mat
load distributional_irfs_4.mat

% Create time index
T = 40;
vTime = linspace(1,T,T);


%%%
% Impulse responses to TFP shock
%%%

% Aggregate variables
figure

subplot(2,3,1)
title('Aggregate Output','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vOutput_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vOutput_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vOutput_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
h = legend('$n_{g} = 2$','$n_{g} = 3$','$n_{g} = 4$');
set(h,'interpreter','latex','fontsize',12)
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,2)
title('Aggregate Consumption','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vConsumption_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vConsumption_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vConsumption_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,3)
title('Aggregate Investment','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vInvestment_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vInvestment_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vInvestment_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,4)
title('Aggregate Hours','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vHours_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vHours_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vHours_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,5)
title('Real Wage','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vWage_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vWage_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vWage_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,6)
title('Real Risk-Free Rate','interpreter','latex','fontsize',14)
hold on
plot(vTime,vInterestRate_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vInterestRate_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vInterestRate_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off


% Distributional variables
figure

subplot(2,2,1)
title('Mean log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_2_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_2_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_2_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
h = legend('$n_{g} = 2$','$n_{g} = 3$','$n_{g} = 4$');
set(h,'interpreter','latex','fontsize',12,'location','northeast')
hold off

subplot(2,2,2)
title('Covariance of productivity w/ log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_4_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_4_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_4_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,2,3)
title('Variance log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_5_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_5_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_5_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,2,4)
title('Marginal utility','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vLog_marginal_utility_TFP_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vLog_marginal_utility_TFP_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vLog_marginal_utility_TFP_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off



%%%
% Investment-Specific Shock
%%%

% Aggregate variables
figure

subplot(2,3,1)
title('Aggregate Output','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vOutput_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vOutput_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vOutput_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
h = legend('$n_{g} = 2$','$n_{g} = 3$','$n_{g} = 4$');
set(h,'interpreter','latex','fontsize',12)
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,2)
title('Aggregate Consumption','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vConsumption_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vConsumption_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vConsumption_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,3)
title('Aggregate Investment','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vInvestment_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vInvestment_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vInvestment_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,4)
title('Aggregate Hours','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vHours_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vHours_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vHours_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,5)
title('Real Wage','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vWage_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vWage_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vWage_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,3,6)
title('Real Risk-Free Rate','interpreter','latex','fontsize',14)
hold on
plot(vTime,vInterestRate_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vInterestRate_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vInterestRate_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off

% Distributional variables
figure

subplot(2,2,1)
title('Mean log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_2_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_2_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_2_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
h = legend('$n_{g} = 2$','$n_{g} = 3$','$n_{g} = 4$');
set(h,'interpreter','latex','fontsize',12,'location','northeast')
hold off

subplot(2,2,2)
title('Covariance of productivity w/ log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_4_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_4_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_4_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,2,3)
title('Variance log capital','interpreter','latex','fontsize',14)
hold on
plot(vTime,vMoment_5_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,vMoment_5_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,vMoment_5_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Deviation','interpreter','latex','fontsize',14)
hold off

subplot(2,2,4)
title('Marginal utility','interpreter','latex','fontsize',14)
hold on
plot(vTime,100*vLog_marginal_utility_Q_Order_2(1:T),'linewidth',1.5,'linestyle','-','color',[178/255,34/255,34/255])
plot(vTime,100*vLog_marginal_utility_Q_Order_3(1:T),'linewidth',1.5,'linestyle','-','color',[10*8/255,62/255,.6*118/255])
plot(vTime,100*vLog_marginal_utility_Q_Order_4(1:T),'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,zeros(size(vTime)),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex','fontsize',14)
ylabel('Percentage deviation','interpreter','latex','fontsize',14)
hold off
