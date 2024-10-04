// Solve for aggregate dynamics using DYNARE
//
// Thomas Winberry, February 1th, 2018

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters.mod"


// Steady state stored in external file (so don't have to re-compute at each step of estimation)
steadyState 											= load('steadyStateForEstimation');

// Value function coefficients
@#for iState in 1 : nState
	parameters valueCoefficientSS_@{iState};
@#endfor

for iState = 1 : @{nState}
	M_.params(@{nCounter} + iState) = steadyState.vCoefficients(iState);
end

@#define nCounter = nCounter + nState

// Capital adjust
@#for iProd in 1 : nProd
	parameters capitalAdjustSS_@{iProd};
@#endfor

for iProd = 1 : @{nProd}
	M_.params(@{nCounter} + iProd) = steadyState.vCapitalAdjust(iProd);
end

@#define nCounter = nCounter + nProd

// Moments
@#for iMoment in 1 : nMeasureCoefficients
	parameters momentSS_@{iMoment};
@#endfor

for iMoment = 1 : @{nMeasureCoefficients}
	M_.params(@{nCounter} + iMoment) = steadyState.vMoments(iMoment);
end

@#define nCounter = nCounter + nMeasureCoefficients

// Measure coefficients
@#for iCoefficient in 1 : nMeasureCoefficients
	parameters measureCoefficientSS_@{iCoefficient};
@#endfor

for iCoefficient = 1 : @{nMeasureCoefficients}
	M_.params(@{nCounter} + iCoefficient) = steadyState.vParameters(iCoefficient);
end

@#define nCounter = nCounter + nMeasureCoefficients

// Prices
parameters wageSS marginalUtilitySS;
M_.params(@{nCounter} + 1) = steadyState.wage;
M_.params(@{nCounter} + 2) = steadyState.marginalUtility;

// Other aggregates
parameters aggregateOutputSS aggregateConsumptionSS aggregateInvestmentSS aggregateHoursSS;
M_.params(@{nCounter} + 3) = steadyState.aggregateOutput;
M_.params(@{nCounter} + 4) = steadyState.aggregateConsumption;
M_.params(@{nCounter} + 5) = steadyState.aggregateInvestment;
M_.params(@{nCounter} + 6) = steadyState.aggregateHours;

// Disutility of labor supply
M_.params(18)               = steadyState.cchi;


//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables.mod"

// Additional variables for estimation
var logAggregateConsumptionObserved logAggregateInvestmentObserved;


//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;
	
	// Model equations
	@#include "equations.mod"
	
	// Additional equations for estimation
	logAggregateConsumptionObserved		= logAggregateConsumption - log(aggregateConsumptionSS);
	logAggregateInvestmentObserved		= logAggregateInvestment - log(aggregateInvestmentSS);

end;


//----------------------------------------------------------------
// Load in steady state externally
// (so don't have to recompute steady state at each step of estimation)
//----------------------------------------------------------------

steady_state_model;

	@#for iState in 1 : nState
		valueCoefficient_@{iState} 				= valueCoefficientSS_@{iState};
	@#endfor
	
	@#for iProd in 1 : nProd
		capitalAdjust_@{iProd}						= capitalAdjustSS_@{iProd};
	@#endfor
	
	@#for iCoefficient in 1 : nMeasureCoefficients
		moment_@{iCoefficient} 					= momentSS_@{iCoefficient};
		measureCoefficient_@{iCoefficient} 	= measureCoefficientSS_@{iCoefficient};
	@#endfor
	
	wage 										= wageSS;
	marginalUtility 							= marginalUtilitySS;
	aggregateOutput 						= aggregateOutputSS;
	aggregateConsumption 				= aggregateConsumptionSS;
	aggregateInvestment 				= aggregateInvestmentSS;
	aggregateHours 						= aggregateHoursSS;
	expectedMarginalUtilityPrime 		= marginalUtilitySS;
	realInterestRate 						= 100 * ((1 / bbeta) - 1);
	
	logAggregateOutput 					= log(aggregateOutputSS);
	logAggregateConsumption 			= log(aggregateConsumptionSS);
	logAggregateInvestment			= log(aggregateInvestmentSS);
	logAggregateHours 					= log(aggregateHoursSS);
	logWage 									= log(wageSS);
	logMarginalUtility						= log(marginalUtilitySS);
	
	aggregateTFP 											= 0;
	aggregateQ 												= 0;
	logAggregateConsumptionObserved 			= 0;
	logAggregateInvestmentObserved 				= 0;
	
end;


//----------------------------------------------------------------
// Set options
//----------------------------------------------------------------

// Specify shock process (=1 to include shock; =0 to not include shock)
shocks;
    var aggregateTFPShock 	= 1;
    var aggregateQShock   	= 1;
end;


// Do not check steady state (un-comment for speed once you know steady state is computed correctly)
options_.steadystate.nocheck = 1;

// Check steady state (un-comment to check whether steady state from .m file satisfies equations.mod)
//steady;
//resid(1)


// Check dynamic regularity conditions (comment out for speed)
//check;
//model_diagnostics;
//model_info;

// Simulate
stoch_simul(order=@{perturbation_order},hp_filter=100,periods=1000,pruning) aggregateTFP logAggregateOutput logAggregateInvestment 
	logAggregateConsumption logAggregateHours logWage realInterestRate moment_2 moment_4 moment_5 logMarginalUtility; 
	
// Estimation
varobs logAggregateConsumptionObserved logAggregateInvestmentObserved;
estimated_params;
rrhoTFP, beta_pdf, 0.9, 0.07;
rrhoQ, beta_pdf, 0.9, 0.07;
ssigmaTFP, inv_gamma_pdf, 0.01, 1;
ssigmaQ, inv_gamma_pdf, 0.01, 1;
corrTFPQ, uniform_pdf, , , -.05, .05;
end;
estimation(order=1,datafile=mDetrendedData,mh_replic=10000,mh_jscale=.8);

// Change directory
cd('../../')
