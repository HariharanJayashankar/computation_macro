% Sets parameter values for replicating figures (steadyStateResults.m)
%
% Thomas Winberry,February 14th, 2018

%----------------------------------------------------------------
% Set model parameters
%----------------------------------------------------------------

global ttheta nnu ddelta rrhoProd ssigmaProd aaUpper aaLower ppsiCapital ...
	bbeta ssigma pphi nSS cchi

% Technology
ttheta 			= .256;								% capital coefficient
nnu 				= .64;								% labor coefficient
ddelta 			= .085;								% depreciation (annual)
rrhoProd 		= .859; 								% persistence of idiosyncratic shocks (annual)
ssigmaProd 	= .022;								% SD innovations of idiosycnratic shocks (annual)
aaUpper 		= .011; 								% no fixed cost region upper bound
aaLower 		= -.011;								% no fixed cost region lower bound
ppsiCapital 	= .0083;							% upper bound on fixed adjustment cost draws

% Preferences
bbeta 			= .961;								% discount factor (quarterly)
ssigma 			= 1;									% coefficient of relative risk aversion
pphi 				= 1 / 1e5;							% inverse Frisch elasticity of labor supply
nSS 				= 1 / 3;								% hours worked in steady state
cchi				= 1;									% labor disutility (will be calibrated in Dynare's steady state file to ensure labor supply = nSS)

% Aggregate shocks (irrelevant for steady state, but included for completeness)
rrhoTFP			= 0.859;							% persistence (annual)
ssigmaTFP		= .014;								% SD of innovations (annual)

%----------------------------------------------------------------
% Set parameters governing approximations
%----------------------------------------------------------------

global nProd nCapital nState prodMin prodMax capitalMin capitalMax nShocks nProdFine nCapitalFine nStateFine ...
	maxIterations tolerance acc dampening nMeasure nStateQuadrature nMeasureCoefficients nProdQuadrature ...
	nCapitalQuadrature kRepSS wRepSS

% Compute representative agent steady state (used in constructing the grids)
global kRepSS wRepSS
kRepSS 			= ((ttheta * (nSS ^ nnu)) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - ttheta));
wRepSS 		= (kRepSS .^ ttheta) * nnu * (nSS ^ (nnu - 1));

% Order of approximation of value function
nProd 			= 3;										% order of polynomials in productivity
nCapital 		= 5;										% order of polynomials in capital
nState 			= nProd * nCapital;					% total number of coefficients

% Bounds on grid space
prodMin 		= -3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
prodMax 		= 3 * ssigmaProd / sqrt(1 - rrhoProd ^ 2);
capitalMin 		= .1 * (exp(prodMin) ^ (1 / (1 - ttheta))) * kRepSS;
capitalMax 		= 2.5 * (exp(prodMax) ^ (1 / (1 - ttheta))) * kRepSS;

% Shocks 
nShocks 		= 3;										% order of Gauss-Hermite quadrature over idiosyncratic shocks

% Finer grid for analyzing policy functions and computing histogram
nProdFine 		= 60;
nCapitalFine 	= 40;
nStateFine 	= nProdFine * nCapitalFine;

% Iteration on value function
maxIterations	= 100;
tolerance 		= 1e-6;
acc 				= 500;									% number of iterations in "Howard improvement step"
dampening 	= 0;							

% Approximation of distribution
nMeasure 					= iMeasureOrder;		% order of polynomial approximating distribution
nProdQuadrature 		= 20; 						% number of quadrature points in productivity dimension
nCapitalQuadrature 	= 20;						% number of quadrature points in capital dimension
nStateQuadrature 		= nProdQuadrature * nCapitalQuadrature;
nMeasureCoefficients 	= (nMeasure * (nMeasure + 1)) / 2 + nMeasure;