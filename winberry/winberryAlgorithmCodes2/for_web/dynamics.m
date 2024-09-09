% Computes and analyzes aggregate dynamics
%
% Thomas Winberry, February 14th, 2018

clear all
close all
clc

cd('./Auxiliary Functions/Dynamics');


%----------------------------------------------------------------
% Set parameters 
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


%----------------------------------------------------------------
% Run Dynare
%----------------------------------------------------------------

eval(['dynare dynamicModel.mod -DnMeasure=' num2str(nMeasure) ' -Dperturbation_order=' num2str(1)])
