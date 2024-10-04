This archive contains MATLAB and DYNARE software for solving the Khan and Thomas (2008) model, using the method
developed in "A Method for Solving and Estimating Heterogeneous Agent Macro Models."  

Thomas Winberry, February 21st, 2018.  Please email me at Thomas.Winberry@chicagobooth.edu with any feedback,
questions, or bug reports.  If you use these codes, please cite: Winberry (2018), "A Method for Solving 
and Estimating Heterogeneous Agent Macro Models."


The following items are included:

--------------------------------------------------------------------------------------------------------
1.  MATLAB software for computing and analyzing the stationary equilibrium of the model with no aggregate
	shocks.  
	

	"steadyState.m"

	This script will solve for the steady state of the model, plot decision rules, plot the stationary distribution,
	compare the parametric family approximation of the distribution to a nonparametric histogram, and compute
	steady state aggregates.  See Appendix A.2 of the paper for a discussion of how to solve the steady state.

	The script calls subroutines:

	a. "set Parameters.m": sets parameters.  This includes both parameters of the economic model, and
		parameters controlling the approximation of the model.
	b.  "coreSteadyState.m": solves for the market-clearing condition in the labor market.  Initially uses a histogram
		approximation of the distribution as an initial guess, and then uses the parametric family approximation of
		the distribution.
	c. "computeGrids.m": creates various grids used in the approximation of the steady state.
		d. "scaleUp.m": transforms [-1,1] to [a,b]; used for computing Chebyshev polynomials
		e. "scaleDown.m": inverse of scaleUp.m
		f. "computeChebyshev.m": compute Chebyshev polynomials over a pre-defined grid.  Used to approximate value function and
			capital accumulation decision conditional on adjusting.
		g. "computeChebyshev2.m": compute Chebyshev polynomials of the 2nd type (derivatives of Chebyshev polynomials) over
			a pre-defined grid.  Used to approximate derivative of the value function.
		h. "computeGaussHermiteQuadrature.m": computes weights and nodes for Gauss-Hermite
			quadrature used to integrate idiosyncratic shocks.  
		i. "computeGaussLegendreQuadrature": computes weights and notes for Gauss-Legendre 
			quadrature used to integrate the parametric family.  
	j. "computePolynomials.m": precomputes Chebyshev polynomials over grids created by "computeGrids.m."  
	k. "computeLMCResidualHistogram.m": computes the difference between specified steady state labor supply nSS and aggregate labor
		demand by firms, using histogram to approximate the distribution.  
		l. "updateCoefficients.m": performs a step of value function iteration.  Given coefficients of the Chebyshev polynomial
			approximating the value function, this function computes optimal policy functions and updates value function	
			"acc" times given those policy functions, where "acc" is a parameter set in "setParameters.m."
		m. "computePolicies.m": given coefficients on capital accumulation conditional on adjusting, computes all investment
			policy rules over a given grid.
		n. "computeDiscreteTransitionMatrix": computes transition matrix defining the evolution of the distribution, given 
			policy rules.  Follows Young (2010).
	o. "computeLMCResidualPolynomials.m:" computes the difference between specified steady state labor supply nSS and aggregate labor
		demand by firms, using parametric family to approximate the distribution.  
		p.  "parametersResidual.m": evaluates the objective function in the minimization problem defining the 
			parameters of the distribution, given moments.

	--------------------------------------


	"steadyStateResults.m"

	This script will replicate Figure 1 and Table 2 from the main text.  The plot is generated automatically.  The entries of Table 2 
	are printed out in vectors.  In addition to the subroutines above, this script calls:

	a. "setParametersFigures.m": sets parameters. This includes both parameters of the economic model, and
		parameters controlling the approximation of the model.  Differs from "setParameters.m" in order to allow looping over
		nMeasureOrder, the order of the polynomial approximating the distribution.



--------------------------------------------------------------------------------------------------------
2.	MATLAB and DYNARE software for computing and analyzing the dynamics of the model with aggregate
	shocks
	
	"dynamics.m"

	This script solves for a local approximation of the model's dynamics using DYNARE.  DYNARE will automatically print
	time-series statistics and plot impulse responses.  Further options can be set in dynamicModel.mod (see below).  In addition to
	subroutines defined above for the steady state, this script calls:
	
	a.  "dynamicModel_steadystate.m": a file to compute the steady state within a DYNARE call; formatted to be in the form 
		required by DYNARE.  Calibrates disutility of labor supply cchi to ensure steady state labor supply = nSS, the
		value specified in "setParameters.m."
	b.  "dynamicModel.mod": the main DYNARE .mod file.  This file performs the following tasks.  First, it defines the parameters
		of the model using "parameters.mod;" importantly, these parameters include the grids and polynomials used to approximate
		the value function, individual decisions, and distribution.  Second, the file defines the variables of the model using 
		"variables.mod."  Third, it defines the model equations in "equations.mod."  Fourth, the file performs the computations 
		necessary to solve the model.  Options are included to check whether the steady state computed by "dynamicModel_steadystate.m"
		actually satisfies the steady state of "equations.mod" and to check that the dynamic regularity conditions necessary
		to compute the solution are satisfied.  It is good practice to check these conditions, but these options can be turned off for speed.
		c.  "parameters.mod": defines the parameters of the model.
		d.  "variables.mod": defines the variables in the model.
		e.  "equations.mod": defines the model equations.

	NB: These .mod files use DYNARE's macro-processor (macro-processor commands are preceeded by @#).  See DYNARE's documentation if you
		are unfamiliar with the macro-processor.  The macro-processor allows the user to write loops in DYNARE.



	--------------------------------------

	"dynamicsResultsFirstOrder.m"

	This script will replicate Figure 2, Table 3, Figure 3, Table 4, and Figure 6 of the paper.  Figures 2, 3, and 6 are generated
	automatically.  The entries of Tables 3 and 4 are printed out in DYNARE output (NB: must turn off investment-specific shocks to 
	replicate the tables).


	--------------------------------------

	"dynamicsResultsSecondOrder.m"

	This script will replicate Figure 4 of the paper.  The figure is generated automatically.  (NB: must turn off investment-specific shocks
	to replicated the figure exactly).


	--------------------------------------

	"dynamicsResultsEstimation.m"

	This script will replicate Table 6, Table 7, and Figure 7 of the paper.  Must set ppsiCapital manually in "dynamicModel_estimation.mod."
	To compute the entries of Table 7, run stoch_simul for estimated parameter values.  In addition to the subroutines described above,
	this script will call:

	a. "dynamicModel_estimation.mod": the main DYNARE .mod file.  It has two key differences from "dynamicModel.mod."  First, it calls the
		estimation routines of DYNARE.  Second, it pre-computes the steady state once and saves it as a parameter, which is then called 
		in the solution using a steady_state_model block.  This greatly speeds up the algorithm since the steady state does not depend 
		on the parameters being estimated.


	--------------------------------------

	"dynamicsResultsAggregation.m"

	This script will replicate Figure 5 and Table 8 of the paper.  The figure is generated automatically.  The table entry corresponds to the 
	case with ssigma_q = vSsigmaQRange(1,1).


