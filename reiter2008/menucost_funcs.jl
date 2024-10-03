using QuantEcon, LinearAlgebra, Roots, Parameters
using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using SparseArrays, SparsityDetection
using Optim

# find stationary distribution
# of markov process characterized by
# transition matrix P
function findStationary(P;tol=1e-6, maxiter=1000)

    nstates = size(P, 1)
    P0 = ones(nstates, 1) ./ nstates
    error=10
    iter = 0

    while iter < maxiter && error > tol
        P1 = P' * P0 

        error = maximum(abs.(P1 - P0))
        iter += 1
        P0 = P1

    end

    if iter == maxiter
        println("Warning. Stationary markov dist: Max iterations reached.")
    end

    return P0

end

#==
Flex price steady state solution
==#
function flexsoln(ϵ, agrid, aPstationary, ζ, ν )

    # get flexi price solution (which is the static solution)
    # w = wage
    # outputs a price vector of same sizee as the shock grid

    function aggPi(L)
        return L * ( (ϵ/(ϵ-1))^(1-ϵ) - (ϵ/(ϵ-1))^(-ϵ)  ) * wout^(1-ϵ) * Pdf^(ϵ-2)
    end
    
    function consumer_lab(L)
        return wout / (wout * L - aggPi(L)) - ζ*L^(1/ν)
    end

    Pdf = (sum(  (exp.(agrid)) .^ (ϵ-1) .* aPstationary))^(1/(ϵ-1)) # price disp
    wout = ( (ϵ-1)/ϵ ) * (1/Pdf)
    pout = sum(((ϵ/(ϵ - 1)) .* wout ./ (exp.(agrid))) .* aPstationary);

    Lout = find_zero(consumer_lab, 1)
    Yout = Lout/Pdf

    return pout, Pdf, wout, Lout, Yout
    
end

#==
value function iterations
v(p, a)
where p is firms last price
a is shock value
Gamma is aggregate state
assume the steady state equilibrium we have is s.t.
aggregate dist doesnt change, therefore all the firm
needs is aggregate Y and wages, and not the whole joint distribution
of p and a
agg is a struct with fields w and Y
==#
function  viterFirm(agg, params; 
                    stochdiscfactor = 1,
                    maxiter=10000, tol=1e-6, 
                    printinterval=1000, printinfo=true)


    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params

    # initial values for v
    vadjust = zeros(np, na)
    vnoadjust = zeros(np, na)
    polp = zeros(Int64, 1, na)
    pollamb = BitArray(undef, np, na)
    pstar = Matrix{CartesianIndex{2}}(undef, 1, na)
    v1 = zeros(np, na)

    # preallocate profit matrix
    profit_mat = zeros(params.np, params.na);
    for pidx=1:np
        for aidx=1:na
            pval = params.pgrid[pidx];
            aval = params.agrid[aidx];
            profit_mat[pidx, aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    error = 10;
    iter = 0;

    while (error > tol) && (iter < maxiter)

        v0 = v1
        ev0 = v0 * aP';

        # iterate over choices
        vnoadjust = profit_mat + β * stochdiscfactor * ev0
        vadjust_val = profit_mat .- κ +  β * stochdiscfactor * ev0
        vadjustmax, pstar = findmax(vadjust_val, dims=1)

        vadjust = repeat(vadjustmax, np, 1)

        v1 = max.(vadjust, vnoadjust)

        error = maximum(abs.(v1 - v0))
        iter += 1;


        if (iter == 1 || mod(iter, printinterval) == 0) && printinfo
            println("Iterations: $iter, Error: $error")
        end
        

    end

    if iter == maxiter
        println("Warning: max iterations reached. Problem may have not converged")
    end

    if printinfo
        println("Firm Viter: Final iterations $iter, Final error $error")
    end

   
    pollamb = vadjust .> vnoadjust
    vout = max.(vadjust, vnoadjust)

    # adjusting policies are jsut one dimensional
    polp = getindex.(pstar, 1)
    polp = polp[:]
    vadjust = vadjust[1, :]

    # interpolate the policies
    itp = interpolate(p.pgrid[polp], (BSpline(Cubic())))
    eitp = extrapolate(itp, Line())
    polp_interp = Interpolations.scale(eitp, p.agrid)

    itp = interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
    eitp = extrapolate(itp, Line())
    pollamb_interp = Interpolations.scale(eitp, p.pgrid, p.agrid)

    return vout, vadjust, vnoadjust, polp_interp, pollamb_interp, iter, error

end

#==
Backward induction of V - for when we know tomorrow's value function 
infl = log(P_t/P_{t-1})
==#
function  vBackwardFirm(agg, params, Z, v1, infl;
                        stochdiscfactor = 1.0,
                        maxiter=10000, tol=1e-6, 
                        printinterval=1000, printinfo=true) 


    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params
    

    # preallocate profit matrix
    profit_infl_mat = zeros(eltype(v1), params.np, params.na);
    profit_mat = zeros(eltype(v1), params.np, params.na);
    for pidx=1:np
        for aidx=1:na
            pval = params.pgrid[pidx];
            pval_adj = pval / (1.0 + infl)
            aval =  log(Z) + params.agrid[aidx];
            profit_mat[pidx, aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
            profit_infl_mat[pidx, aidx] = (pval_adj^(1-ϵ) - pval_adj^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    error = 10;
    iter = 0;

    # interpolate v1
    itp = interpolate(v1, (BSpline(Cubic()), BSpline(Cubic())))
    eitp = extrapolate(itp, Line())
    v1interp = Interpolations.scale(eitp, pgrid, agrid)
    v1_infl = zeros(eltype(v1), np, na)
    v1_noinfl = zeros(eltype(v1), np, na)
    for pidx = 1:np
        for aidx = 1:na
            pval = pgrid[pidx]
            pval_adj = pval / (1.0 + infl)

            aval = agrid[aidx]
            v1_infl[pidx, aidx] = v1interp(pval_adj, aval)
            v1_noinfl[pidx, aidx] = v1interp(pval, aval)
        end
    end

    ev1_infl = v1_infl * aP';
    ev1_noinfl = v1_noinfl * aP';

    # iterate over choices
    vnoadjust = profit_infl_mat + β * stochdiscfactor * ev1_infl
    vadjust_val = profit_mat .- κ +  β * stochdiscfactor * ev1_noinfl
    vadjustmax, pstar = findmax(vadjust_val, dims=1)
    vadjust = repeat(vadjustmax, np, 1)

    pollamb = vadjust .> vnoadjust
    vout = max.(vadjust, vnoadjust)

    polp = getindex.(pstar, 1)
    polp = polp[:]
    vadjust = vadjust[1, :]

    # interpolate the policies
    itp = interpolate(p.pgrid[polp], (BSpline(Cubic())))
    eitp = extrapolate(itp, Line())
    polp_interp = Interpolations.scale(eitp, p.agrid)

    itp = interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
    eitp = extrapolate(itp, Line())
    pollamb_interp = Interpolations.scale(eitp, p.pgrid, p.agrid)

    return vout, vadjust, vnoadjust, p.pgrid[polp], pollamb_interp, iter, error

end

#==
Make transition function for distribution of (p,a)
Given the policy functions for price changes
==#
function Tfunc(omega0, polp, pollamb, params)

    aP = params.aP

    # update shock dist
    omega1hat = zeros(eltype(omega0), params.np, params.na);
    for pidx = 1:params.np
        for aidx = 1:params.na

            for a0idx = 1:params.na
                omega1hat[pidx, aidx] = omega1hat[pidx, aidx] + omega0[pidx, a0idx] * aP[a0idx, aidx];
            end
        end
    end

    # update policies
    omega1 = zeros(eltype(omega0), params.np, params.na);
    for pidx = 1:params.np
        for aidx = 1:params.na
            
            p1idx = polp[aidx];

            # non adjusters
            omega1[pidx, aidx] = omega1[pidx, aidx] + (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
            
            # adjusters
            omega1[p1idx, aidx] = omega1[p1idx, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];
        end
    end

    # omega1 = sparse(omega1)
    # omega1hat  = sparse(omega1hat)
    # droptol!(omega1, 1e-10)
    # droptol!(omega1hat, 1e-10)
    return omega1, omega1hat

end

#==
yougn simulation of updating exogenous shock process
==#
function Tfunc_updateshocks(omega0, params, ngrid, Z)

    aP = params.aP

    omega1hat = zeros(ngrid, params.na);
    agrid = log(Z) .+ params.agrid
    for pidx = 1:ngrid
        for aidx = 1:params.na

            aval = agrid[aidx]

            if aval > params.agrid[1] && aval < params.agrid[end]
                aidx_vals = searchsorted(params.agrid, aval)
                aidx_lo = last(aidx_vals)
                aidx_hi = aidx_lo + 1
                total_dist = params.agrid[aidx_hi] - params.agrid[aidx_lo]

                wt_lo = 1.0 - (aval - params.agrid[aidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo

                for a0idx = 1:params.na
                    omega1hat[pidx, aidx_hi] = omega1hat[pidx, aidx_hi] + wt_hi * omega0[pidx, a0idx] * aP[a0idx, aidx_hi]
                    omega1hat[pidx, aidx_lo] = omega1hat[pidx, aidx_lo] + wt_lo * omega0[pidx, a0idx] * aP[a0idx, aidx_lo]
                end

                

            elseif aval <= params.agrid[1]

                for a0idx = 1:params.na
                    omega1hat[pidx, 1] = omega1hat[pidx, 1] +  omega0[pidx, a0idx] * aP[a0idx, 1]
                end

            elseif aval >= params.agrid[end]

                for a0idx = 1:params.na
                    omega1hat[pidx, end] = omega1hat[pidx, end] +  omega0[pidx, a0idx] * aP[a0idx, end]
                end

            end
        end
    end

    return omega1hat

end

#==
Tfunc general is meant to be a typical Young 2010 non stochastic simulation
function.
- omega0 is a 2d matrix of distribvution over (p,a) in the previous period
- omega0hat is 2d matrix of distribution over (p,a) at beginin of period
after shocks have been realized but before pricing decisions have been made
- polp is the policy function conditional oin adjustment - crucially in this
function it returns the actual price value instead of the price index in the
pricegrid
- pollamb is {1,0} for whether a firm adjusts or not

If inflation is non zero you can apss in the discretized policy funciton matrices
which spit out the index of the pgrid

If inflation is non zero use should pass in the interpolated policy functions
which take in the actual pvalues and spit out the pvalues aswell
==#
function Tfunc_general(omega0, polp, pollamb, params, ngrid, infl, Z)
    

    omega1hat = Tfunc_updateshocks(omega0, params, ngrid, Z)

    # update policies
    omega1 = zeros(ngrid, params.na);
    

    for pidx = 1:ngrid
        for aidx = 1:params.na

            # non adjusters
            pval0 = params.pgrid_fine[pidx]
            pval = pval0 / (1.0 + infl)
            aval = log(Z) + params.agrid[aidx]
            if pval0 > params.plo && pval0 < params.phi
                pidx_vals = searchsorted(params.pgrid_fine, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid_fine[pidx_hi] - params.pgrid_fine[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid_fine[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1 - wt_lo

                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * (1.0 - pollamb(pval0, aval)) * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * (1.0 - pollamb(pval0, aval)) * omega1hat[pidx, aidx];
                

            elseif pval <= params.plo

                omega1[1, aidx] = omega1[1, aidx] +  (1.0 - pollamb(pval0, aval)) * omega1hat[pidx, aidx];

            elseif pval >= params.phi

                omega1[end, aidx] = omega1[end, aidx] +  (1.0 - pollamb(pval0, aval)) * omega1hat[pidx, aidx];

            end

            # adjusters
            pval0 = params.pgrid_fine[pidx]
            pval = polp(aval)
            if pval > params.plo && pval < params.phi
                pidx_vals = searchsorted(params.pgrid_fine, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid_fine[pidx_hi] - params.pgrid_fine[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid_fine[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo
                
                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * pollamb(pval0, aval) * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * pollamb(pval0, aval) * omega1hat[pidx, aidx];

            elseif pval <= params.plo

                omega1[1, aidx] = omega1[1, aidx] + pollamb(pval0, aval) * omega1hat[pidx, aidx];

            elseif pval >= params.phi

                omega1[end, aidx] = omega1[end, aidx] + pollamb(pval0, aval) * omega1hat[pidx, aidx];

            end
        end
    end

    # omega1 = omega1 ./ sum(omega1)
    # omega1hat = omega1hat ./ sum(omega1hat)
    return omega1, omega1hat

end


function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)


    omega1 = ones(params.np_fine, params.na);
    omega1 = omega1 ./ (params.np_fine*params.na);
    omega1hat = zeros(params.np_fine, params.na);
    # omega1 = sparse(omega1)
    # omega1hat = sparse(omega1hat)
    error = 10;
    iter = 0;
    
    while (error > tol) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc_general(omega0, polp, pollamb, params, params.np_fine, params.Π_star, 1.0)
        error = maximum(abs.(omega1 - omega0))
        iter += 1;
        omega0hat = omega1hat;

        if ((iter==1) || (mod(iter, printinterval) == 0)) & printinfo
            println("Iterations: $iter, Error $error")
        end

    end

    if printinfo
        println("Final Iterations: $iter, Final error: $error")
    end

    return omega1hat, omega1

end

#==
Find equilibiurm Y and w to clear steady state
==#
function findEquilibrium_ss(p; winit=1, tol=1e-4, max_iter=200, deltaw=0.1,
                            Yinit=1, deltaY=0.1,
                            printinterval=10)

    # gettngg value funcion
    w0 = winit;
    Y0 = Yinit;
    iter = 0;
    error = 10

    # preallocating
    Vadjust = zeros(p.na)
    Vnoadjust = zeros(p.np, p.na)
    polp = zeros(p.na)
    pollamb = zeros(p.np, p.na)
    omega = zeros(p.np, p.na)
    omegahat = zeros(p.np, p.na)
    C = 0.0

    # outer labour market loop
    while iter < max_iter && error > tol

        agg = (w=w0, Y=Y0);
        V, Vadjust ,Vnoadjust, polp, pollamb  = viterFirm(agg, p; maxiter=10000, tol=1e-6, printinfo=false)

        # get joint distribution of prices and shocks
        omegahat, omega = genJointDist(polp, pollamb, p; printinfo=false);

        # get implied aggregate price
        aggprice = 0.0
        for pidx=1:p.np
            for aidx = 1:p.na
                pval = p.pgrid[pidx]
                pchange = pollamb[pidx, aidx]
                p1val = pchange * p.pgrid[polp[aidx]]  + (1.0 - pchange) * pval
                aggprice += p1val^(1.0 - p.ϵ) * omegahat[pidx, aidx]
            end
        end
        aggprice = aggprice^(1.0/(1.0-p.ϵ))

        # get profits to give HH
        # get aggregate fixed cost payments
        # get labour demand
        # integrate
        Ld = 0.0
        F = 0.0
        for pidx = 1:p.np
            for aidx = 1:p.na
                pval = p.pgrid[pidx]
                a1val = p.agrid[aidx]
                pchange = pollamb[pidx, aidx]
                p1val = pchange * p.pgrid[polp[aidx]] + (1.0 - pchange)*pval
                F += p.κ * pchange * omegahat[pidx, aidx]
                Ld += p1val^(-p.ϵ) * exp(-a1val) * Y0 * omegahat[pidx,aidx]
            end
        end

        C = Y0 - F
        w_implied = p.ζ * Ld^(1/p.ν) * C

        # updating guesses
        errorY = abs(aggprice - 1.0)
        errorw = abs(w0 - w_implied)
        error = max(errorY, errorw)

        if aggprice > 1.0
            Y1 = Y0 * (1.0 - deltaY)
        else
            Y1 = Y0 * (1.0 + deltaY)
        end
        Y1 = (1-deltaY) * Y0 + deltaY * Yimplied;
        w1 = (1-deltaw) * w0 + deltaw * w_implied;
        Y0 = Y1
        w0 = w1

        iter += 1

        if iter == 1 || mod(iter, printinterval) == 0
            println("Iterations: $iter")
            println("Error w: $errorw, W guess: $w0")
            println("Error Y: $errorY, Y guess: $Y0")
        end

    end

    return w0, Y0, Vadjust, Vnoadjust, polp, pollamb,  omega, omegahat, C, iter, error
end

#==
x = [w , y]
==#
function equilibriumResidual(x, p)

    w = x[1]
    Y = x[2]
    agg = (w=w, Y=Y);
    V, Vadjust ,Vnoadjust, polp, pollamb  = viterFirm(agg, p; maxiter=10000, tol=1e-6, printinfo=false)
    itp = interpolate(polp, (BSpline(Cubic())))
    eitp = extrapolate(itp, Line())
    polp_l_interp = Interpolations.scale(eitp, p.agrid)

    # get joint distribution of prices and shocks
    omegahat, omega = genJointDist(polp, pollamb, p; printinfo=false);

    # get implied aggregate price
    aggprice = 0.0
    for pidx=1:p.np_fine
        for aidx = 1:p.na
            pval = p.pgrid_fine[pidx]
            aval = p.agrid[aidx]
            pchange = pollamb(pval, aval)
            p1val = pchange * polp(aval)  + (1.0 - pchange) * pval
            aggprice += p1val^(1.0 - p.ϵ) * omegahat[pidx, aidx]
        end
    end
    aggprice = aggprice^(1.0/(1.0-p.ϵ))

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0.0
    F = 0.0
    for pidx = 1:p.np_fine
        for aidx = 1:p.na
            pval = p.pgrid_fine[pidx]
            a1val = p.agrid[aidx]
            pchange = pollamb(pval, a1val)
            p1val = pchange * polp(a1val) + (1.0 - pchange)*pval
            F += p.κ * pchange * omegahat[pidx, aidx]
            Ld += p1val^(-p.ϵ) * exp(-a1val) * Y * omegahat[pidx,aidx]
        end
    end

    C = Y - F
    w_implied = p.ζ * Ld^(1/p.ν) * C

    # sum of squared errors
    errorY = (aggprice - 1.0)^2.0
    errorw = (w - w_implied)^2.0
    error = errorY + errorw
    return error, w, Y, Vadjust, Vnoadjust, polp, pollamb, omega, omegahat, C

end




# Solving using Reiter 2009

#==
Residual equations which equal 0 at equilibrium
Xl is lagged value of variables
X is current value of variables
The variables are:
- each bin of (p,a) to track distribution - stacked as a vector, size np * na
- w
- r
- Y
- C
- Z (aggregate shock)
- Ey
- Ew
- EMu(C) - expected marginal utility
- V - value functions

A lot of it is rewriting the findEquilibrium_ss function
but without iterating till steady state for the distribution
Note r and C are added in compared to steady state since
the euler equation has to hold in equilibrium outside
SS now
Last three are needed to solve todays value function

==#
function residequations(Xl, X, 
                η, ϵ, 
                p, yss )
    
    # default params
    @unpack np, na, np_fine = p

    sizedist = np_fine * na
    sizev = np * na

    # unpacking
    omega_l = reshape(Xl[1:sizedist], np_fine, na)
    omega = reshape(X[1:sizedist], np_fine, na)

    Vl = reshape(Xl[(sizedist+1):(sizedist+sizev)], np, na)
    V = reshape(X[(sizedist+1):(sizedist+sizev)], np, na)

    wl, rl, Yl, Cl, Zl  = Xl[(sizedist+sizev+1):(end-1)]
    w, r, Y, C, Z = X[(sizedist+sizev+1):(end-1)]
    infl_l = Xl[end]
    infl = X[end]

    # expectation errors
    ηv = reshape(η[1:sizev], np, na)
    # η_polp = η[(na+1):2*na]
    # ηv_noadj = reshape(η[(2*na+1):(2*na + sizev)], np, na)
    # ηv_noadj = reshape(η[(na+1):(na + sizedist)], np, na)
    η_ee = η[sizev + 1]


    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    stochdiscfactor = Cl/C
    # calculate implied polp
    V_l_check, Vadj_l_check, Vnoadj_l_check, polp_l_check, pollamb_l_check, _, _ = vBackwardFirm(
        (Y=Y, w=w), p, Z, V, infl, stochdiscfactor = stochdiscfactor
    )

    V_l_check += ηv
    itp = interpolate(polp_l_check, (BSpline(Cubic())))
    eitp = extrapolate(itp, Line())
    polp_l_interp = Interpolations.scale(eitp, p.agrid)


    #==
    Compute Distribution checks
    ==#
    # this gives distribution at the start for period t before period t shocks
    # have been realized
    omega1, omega0hat = Tfunc_general(omega_l, polp_l_interp, pollamb_l_check, p, p.np_fine, infl, Z)
    # omega1hat = Tfunc_updateshocks(omega1, p, p.np_fine, Z)

    # get implied aggregate Y
    pdist_l = sum(omega1, dims=2)
    aggprice_l = (sum((p.pgrid_fine).^(1.0-p.ϵ) .* pdist_l))^(1.0/(1.0 - p.ϵ))

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0.0
    F = 0.0
    for pidx = 1:p.np_fine
        for aidx = 1:p.na
            pval = p.pgrid_fine[pidx]
            aval = p.agrid[aidx]
            pchange = pollamb_l_check(pval, aval)

            F += p.κ * pchange * omega0hat[pidx, aidx]
            Ld += pval^(-p.ϵ) * exp(-aval) * Y * omega1[pidx,aidx]
        end
    end

    # monetary policy
    # talor rule
    r_val = rl +  p.ϕ_infl*(infl - p.Π_star) + p.ϕ_output*(Y-yss) + ϵ[2]
    mon_pol_error = r_val - r

    cerror = C - Y + F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = 1.0/Cl - (1+rl)*p.β/(C*(1.0+infl)) - η_ee
    zerror = log(Z) - p.ρ_agg * log(Zl) - ϵ[1]


    # == residual vector == #
    residvector = zero(X)
    residvector[1:sizedist] = vec(omega1 - omega)
    residvector[(sizedist+1):(sizedist+sizev)] = vec(V_l_check - Vl)


    # other error
    residvector[(sizedist+sizev+1):end] = [w-w_implied,1.0 - aggprice_l,
                                        euler_error,cerror,
                                        mon_pol_error,
                                        zerror]

     
    return residvector
    
end

# make wrapper functions of only each argument of resid to alloww forward diff
# not get worried about types


#==
Linarizeed coefficients
==#
function linearized_coeffs(equations, xss, ϵ_ss, p)

    # l for lag
    function residequations_lti_l(Xl, xss, ϵ, p)

        xss = convert.(eltype(Xl), xss)
        ϵ = convert.(eltype(Xl), ϵ)
        residout = equations(Xl, xss, xss, ϵ, p)

    end

    # c for current
    function residequations_lti_c(X, xss, ϵ, p)

        xss = convert.(eltype(X), xss)
        ϵ = convert.(eltype(X), ϵ)
        residout = equations(xss, X, xss, ϵ, p)

    end

    # f for future
    function residequations_lti_f(Xf, xss, ϵ, p)

        xss = convert.(eltype(Xf), xss)
        ϵ = convert.(eltype(Xf), ϵ)
        residout = equations(xss, xss, Xf, ϵ, p)

    end

    # eps for epsilon
    function residequations_lti_eps(xss, ϵ, p)

        xss = convert.(eltype(ϵ), xss)
        residout = equations(xss, xss, xss, ϵ, p)

    end


    H1 = ForwardDiff.jacobian(t -> residequations_lti_l(t, xss, ϵ_ss, p), xss)
    H2 = ForwardDiff.jacobian(t -> residequations_lti_c(t,  xss, ϵ_ss, p), xss)
    H3 = ForwardDiff.jacobian(t -> residequations_lti_f(t, xss, ϵ_ss, p), xss)
    H4 = ForwardDiff.jacobian(t -> residequations_lti_eps(xss, t, p), ϵ_ss)

    return H1, H2, H3, H4

end