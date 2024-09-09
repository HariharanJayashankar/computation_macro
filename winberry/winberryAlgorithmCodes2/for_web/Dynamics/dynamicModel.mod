// Solve for aggregate dynamics using DYNARE
//
// Thomas Winberry, February 14th, 2018

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters.mod"


//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables.mod"


//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

	@#include "equations.mod"

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

// Change directory
cd('../../')
