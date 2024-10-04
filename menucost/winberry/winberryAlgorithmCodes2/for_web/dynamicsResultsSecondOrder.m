% Computes and analyzes first-order approximation of model dynamics
%
% Thomas Winberry, February 14th, 2018

cd('./Auxiliary Functions/Dynamics');

%----------------------------------------------------------------
% Set economic parameters 
%----------------------------------------------------------------

setParameters;

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

% Grids
computeGrids;

% Polynomials over grids
computePolynomials;

%----------------------------------------------------------------
% Save parameters for Dynare 
%----------------------------------------------------------------

% Economic parameters
save economicParameters.mat ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS rrhoTFP ssigmaTFP rrhoQ ssigmaQ corrTFPQ cchi
	
% Approximation parameters
save approximationParameters.mat nProd nCapital nState nProdQuadrature nCapitalQuadrature nStateQuadrature ...
	nMeasureCoefficients nMeasure prodMin prodMax capitalMin capitalMax nShocks
	
% Grids
save grids.mat vShocksGrid vShocksWeights mStateGrid mQuadratureGrid vQuadratureWeights mProdPrimeQuadrature

% Polynomials
save polynomials.mat mStatePoly mStatePoly2 vStatePolySquared aProdPrimePoly mQuadraturePoly aProdPrimeQuadraturePoly

% Set gridsizes manually (so that you can pick out relevant series below)
nProd                   		= 3;
nCapital                	= 5;
nShocks                 	= 3;
nState                  		= nProd * nCapital;
nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure;

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

eval(['dynare dynamicModel.mod noclearall -DnMeasure=' num2str(nMeasure) ' -Dperturbation_order=' num2str(2)])

%----------------------------------------------------------------
% Check if impulse response is sign-dependent (follows Pfiefer's GIRF codes)
%----------------------------------------------------------------

% Initialize
T = 40; tBurn = 200; nReplication = 1000;
aIRFPos = zeros(M_.endo_nbr,tBurn+T,nReplication);
aIRFNeg = zeros(M_.endo_nbr,tBurn+T,nReplication);
tic;
% Do the replications
for iReplication = 1 : nReplication
	
	% Shock sequence for current run
	vShocks = randn(tBurn+T,2);
	
	% Positive shock
	vShocksImpulse = vShocks; vShocksImpulse(tBurn+1,1) = vShocksImpulse(tBurn+1,1) + 1;
	vYSeries = simult_(oo_.dr.ys,oo_.dr,vShocks,options_.order);
	vYSeriesImpulse = simult_(oo_.dr.ys,oo_.dr,vShocksImpulse,options_.order);
	aIRFPos(:,:,iReplication) = vYSeriesImpulse(:,2:tBurn+T+1) - vYSeries(:,2:tBurn+T+1);
	
	% Negative shock
	vShocksImpulse = vShocks; vShocksImpulse(tBurn+1,1) = vShocksImpulse(tBurn+1,1) - 1;
	vYSeries = simult_(oo_.dr.ys,oo_.dr,vShocks,options_.order);
	vYSeriesImpulse = simult_(oo_.dr.ys,oo_.dr,vShocksImpulse,options_.order);
	aIRFNeg(:,:,iReplication) = vYSeriesImpulse(:,2:tBurn+T+1) - vYSeries(:,2:tBurn+T+1);
	
end
toc;
mIRFPos = mean(aIRFPos,3); mIRFNeg = mean(aIRFNeg,3);

% Extract relevant series 
vLogOutputPos = mIRFPos(nState+nProd+2*nMeasureCoefficients+9,tBurn+1:tBurn+T)'; vLogOutputNeg = mIRFNeg(nState+nProd+2*nMeasureCoefficients+9,tBurn+1:tBurn+T)';
vLogInvestmentPos = mIRFPos(nState+nProd+2*nMeasureCoefficients+11,tBurn+1:tBurn+T)'; vLogInvestmentNeg = mIRFNeg(nState+nProd+2*nMeasureCoefficients+11,tBurn+1:tBurn+T)';


%----------------------------------------------------------------
% Check if impulse response is state-dependent (follows Pfiefer's GIRF codes)
%----------------------------------------------------------------

% Initialize
T = 40; tBurn = 200; nReplication = 1000;
aIRFGood = zeros(M_.endo_nbr,tBurn+T,nReplication);
aIRFBad = zeros(M_.endo_nbr,tBurn+T,nReplication);
tic;
% Do the replications
for iReplication = 1 : nReplication
	
	% Shock sequence for current run
	vShocksBaseline = randn(tBurn+T,2);
	
	% Good history
	vShocks = vShocksBaseline; vShocks(tBurn-1:tBurn,1) = vShocks(tBurn-1:tBurn,1) + 1;
	vShocksImpulse = vShocks; vShocksImpulse(tBurn+1,1) = vShocksImpulse(tBurn+1,1) + 1;
	vYSeries = simult_(oo_.dr.ys,oo_.dr,vShocks,options_.order);
	vYSeriesImpulse = simult_(oo_.dr.ys,oo_.dr,vShocksImpulse,options_.order);
	aIRFGood(:,:,iReplication) = vYSeriesImpulse(:,2:tBurn+T+1) - vYSeries(:,2:tBurn+T+1);
	
	% Negative shock
	vShocks = vShocksBaseline; vShocks(tBurn-1:tBurn,1) = vShocks(tBurn-1:tBurn,1) - 1;
	vShocksImpulse = vShocks; vShocksImpulse(tBurn+1,1) = vShocksImpulse(tBurn+1,1) + 1;
	vYSeries = simult_(oo_.dr.ys,oo_.dr,vShocks,options_.order);
	vYSeriesImpulse = simult_(oo_.dr.ys,oo_.dr,vShocksImpulse,options_.order);
	aIRFBad(:,:,iReplication) = vYSeriesImpulse(:,2:tBurn+T+1) - vYSeries(:,2:tBurn+T+1);
	
end
toc;
mIRFGood = mean(aIRFGood,3); mIRFBad = mean(aIRFBad,3);

% Extract relevant series
vLogOutputGood = mIRFGood(nState+nProd+2*nMeasureCoefficients+9,tBurn+1:tBurn+T)'; vLogOutputBad = mIRFBad(nState+nProd+2*nMeasureCoefficients+9,tBurn+1:tBurn+T)';
vLogInvestmentGood = mIRFGood(nState+nProd+2*nMeasureCoefficients+11,tBurn+1:tBurn+T)'; vLogInvestmentBad = mIRFBad(nState+nProd+2*nMeasureCoefficients+11,tBurn+1:tBurn+T)';


%----------------------------------------------------------------
% Check if impulse response is state-dependent (follows Pfiefer's GIRF codes)
%----------------------------------------------------------------

% Plot
vTime = linspace(1,T,T);

figure

subplot(2,2,1)
hold on
plot(vTime,100*vLogOutputPos,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,-100*vLogOutputNeg,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
h = legend('Positive shock','Negative shock');
set(h,'interpreter','latex','fontsize',12)
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex')
ylabel('Percentage deviation','interpreter','latex')
title('Sign Dependence, Output','interpreter','latex','fontsize',14)
hold off

subplot(2,2,2)
hold on
plot(vTime,100*vLogInvestmentPos,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,-100*vLogInvestmentNeg,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex')
ylabel('Percentage deviation','interpreter','latex')
title('Sign Dependence, Investment','interpreter','latex','fontsize',14)
hold off

subplot(2,2,3)
hold on
plot(vTime,100*vLogOutputGood,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,100*vLogOutputBad,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
h = legend('Expansion','Recession');
set(h,'interpreter','latex','fontsize',12)
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex')
ylabel('Percentage deviation','interpreter','latex')
title('State Dependence, Output','interpreter','latex','fontsize',14)
hold off

subplot(2,2,4)
hold on
plot(vTime,100*vLogInvestmentGood,'linewidth',1.5,'linestyle','-','color',[8/255,62/255,118/255])
plot(vTime,100*vLogInvestmentBad,'linewidth',1.5,'linestyle','--','color',[178/255,34/255,34/255])
plot(vTime,zeros(T,1),'linewidth',1.5,'linestyle','--','color','k')
set(gcf,'color','w')
grid on
xlim([1 T])
xlabel('Years since shock','interpreter','latex')
ylabel('Percentage deviation','interpreter','latex')
title('State Dependence, Investment','interpreter','latex','fontsize',14)
hold off
