% Estimates model
%
% Thomas Winberry, February 19th, 2018


cd('./Auxiliary Functions/Dynamics');


%----------------------------------------------------------------
% Set up grids
%----------------------------------------------------------------

% Set parameters
setParameters_estimation;

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


%----------------------------------------------------------------
% Pre-compute steady state
% (so don't have to re-compute it at each iteration)
%----------------------------------------------------------------

% Solve for steady state wage
coreSteadyState;

% Compute steady state objects
[vCoefficientsNew,vCapitalAdjust] 		= updateCoefficients(vCoefficients);
[resid,vMoments,vParameters,aggregateConsumption,marginalUtility,aggregateOutput,aggregateCapital,aggregateInvestment] = ...
	computeLMCResidualPolynomials(wage,vParameters,vMomentsHistogram,mGridMoments);
aggregateHours 								= resid + nSS;
cchi 													= (wage * marginalUtility / (aggregateHours ^ pphi));
	
% Rename for compatability with dynare file
vCoefficients										= marginalUtility * vCoefficientsNew;
vParameters										= vParameters(2:nMeasureCoefficients+1);

% Save steady state in .mat file
save steadyStateForEstimation.mat vCoefficients vCapitalAdjust vMoments vParameters wage marginalUtility aggregateConsumption ...
	aggregateOutput aggregateCapital aggregateInvestment aggregateHours cchi

%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

eval(['dynare dynamicModel_estimation.mod noclearall -DnMeasure=' num2str(nMeasure) ' -Dperturbation_order=' num2str(1)])
	
