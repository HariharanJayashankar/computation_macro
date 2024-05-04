using QuantEcon, LinearAlgebra, Roots, Parameters
using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using SparseArrays, SparsityDetection

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

    polp = getindex.(pstar, 1)
    polp = repeat(polp, np, 1)

    return vout, vadjust, vnoadjust, polp, pollamb, iter, error

end

#==
Backward induction of V - for when we know tomorrow's value function 
infl = log(P_t/P_{t-1})
==#
function  vBackwardFirm(agg, params, Z, v1, infl, infl_f; 
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
            pval_adj = pval / exp(infl)
            aval =  log(Z) .+ params.agrid[aidx];
            profit_mat[pidx, aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
            profit_infl_mat[pidx, aidx] = (pval_adj^(1-ϵ) - pval_adj^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    error = 10;
    iter = 0;

    # interpolate v1
    itp = interpolate(v1, (BSpline(Linear()), NoInterp()))
    eitp = extrapolate(itp, Line())
    v1interp = Interpolations.scale(eitp, pgrid, 1:na)
    v1_infl = zeros(eltype(v1), np, na)
    v1_noadj = zeros(eltype(v1), np, na)
    for pidx = 1:np
        for aidx = 1:na
            pval = pgrid[pidx]
            pval_adj = pval / exp(infl_f)
            pval_noadj = pval_adj / exp(infl)
            v1_infl[pidx, aidx] = v1interp(pval_adj, aidx)
            v1_noadj[pidx, aidx] = v1interp(pval_noadj, aidx)
        end
    end

    ev1_infl = v1_infl * aP';
    ev1_noadj_infl = v1_noadj * aP';

    # iterate over choices
    vnoadjust = profit_mat + β * stochdiscfactor * ev1_infl
    vadjust_val = profit_mat .- κ +  β * stochdiscfactor * ev1_infl
    vadjustmax, pstar = findmax(vadjust_val, dims=1)

    vadjust = repeat(vadjustmax, np, 1)

    pollamb = vadjust .> vnoadjust
    vout = max.(vadjust, vnoadjust)

    polp = getindex.(pstar, 1)
    polp = repeat(polp, np, 1)

    return vout, vadjust, vnoadjust, polp, pollamb, iter, error

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
            
            p1idx = polp[pidx, aidx];

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
Tfunc general is meant to be a typical Young 2010 non stochastic simulation
function.
- omega0 is a 2d matrix of distribvution over (p,a) in the previous period
- polp is the policy function conditional oin adjustment - crucially in this
function it returns the actual price value instead of the price index in the
pricegrid
- pollamb is {1,0} for whether a firm adjusts or not
==#
function Tfunc_general(omega0, polp, pollamb, params, ngrid, infl, Z)
    
    aP = params.aP

    # update shock dist
    omega1hat = zeros(params.np, params.na);
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
                wt_hi = 1 - wt_lo

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

    # update policies
    omega1 = zeros(ngrid, params.na);
    for pidx = 1:ngrid
        for aidx = 1:params.na

            # non adjusters
            pval0 = params.pgrid[pidx]
            pval0 = pval0 / exp(infl)
            if pval0 > params.pgrid[1] && pval0 < params.pgrid[end]
                pidx_vals = searchsorted(params.pgrid, pval0)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid[pidx_hi] - params.pgrid[pidx_lo]

                wt_lo = 1.0 - (pval0 - params.pgrid[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1 - wt_lo

                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                

            elseif pval0 <= params.pgrid[1]

                omega1[1, aidx] = omega1[1, aidx] +  (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            elseif pval0 >= params.pgrid[end]

                omega1[end, aidx] = omega1[end, aidx] +  (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            end

            # adjusters
            pval = params.pgrid[polp[pidx, aidx]];
            pval = pval / exp(infl);
            if pval > params.pgrid[1] && pval < params.pgrid[end]
                pidx_vals = searchsorted(params.pgrid, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid[pidx_hi] - params.pgrid[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1 - wt_lo
                
                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * pollamb[pidx, aidx] * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval <= params.pgrid[1]

                omega1[1, aidx] = omega1[1, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval >= params.pgrid[end]

                omega1[end, aidx] = omega1[end, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            end
        end
    end

    omega1 = omega1 ./ sum(omega1)
    omega1hat = omega1hat ./ sum(omega1hat)
    return omega1, omega1hat

end

function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)


    omega1 = ones(params.np, params.na);
    omega1 = omega1 ./ (params.np*params.na);
    omega1hat = zeros(params.np, params.na);
    # omega1 = sparse(omega1)
    # omega1hat = sparse(omega1hat)
    error = 10;
    iter = 0;
    
    while (error > tol) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc_general(omega0, polp, pollamb, params, params.np, params.Π_star, 1.0)
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
    Vadjust = zeros(p.np, p.na)
    Vnoadjust = zeros(p.np, p.na)
    polp = zeros(p.np, p.na)
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
        pdist = sum(omega, dims=2)

        # get implied aggregate Y
        Yimplied = sum((p.pgrid .^ (-p.ϵ))' * pdist)

        # get profits to give HH
        # get aggregate fixed cost payments
        # get labour demand
        # integrate
        Ld = 0.0
        F = 0.0
        for pidx = 1:p.np
            for aidx = 1:p.na
                pval = p.pgrid[pidx]
                aval = p.agrid[aidx]
                F += p.κ * pollamb[pidx, aidx] .* omegahat[pidx, aidx] # who adjusts in a period
                Ld += pval^(-p.ϵ) * exp(-aval) * Y0 * omega[pidx,aidx]
            end
        end

        C = Yimplied - F
        w_implied = p.ζ * Ld^(1/p.ν) * C

        # updating guesses
        errorY = abs(Yimplied - Y0)
        errorw = abs(w0 - w_implied)
        error = max(errorY, errorw)
        # error = errorY

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
    @unpack np, na = p

    sizedist = np * na

    # unpacking
    omega_l = reshape(Xl[1:sizedist], np, na)
    omega = reshape(X[1:sizedist], np, na)

    V_l = reshape(Xl[(sizedist+1):(2*sizedist)], np, na)
    V = reshape(X[(sizedist+1):(2*sizedist)], np, na)

    wl, rl, Yl, Cl, Zl  = Xl[(2*sizedist+1):(end-1)]
    w, r, Y, C, Z = X[(2*sizedist+1):(end-1)]
    infl_l = Xl[end]
    infl = X[end]

    # expectation errors
    ηv = reshape(η[1:sizedist], np, na)
    η_ee = η[end]

    stochdiscfactor = (Cl*exp(infl))/C 

    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    # calculate implied polp
    V_l_check, Vadj_l, Vnoadj_l, polp_l_check, pollamb_l, _, _ = vBackwardFirm(
        (Y=Yl, w=wl), p, Zl, V, infl_l, infl, stochdiscfactor = stochdiscfactor
    )
    V_l_check += ηv

    #==
    Compute Distribution checks
    ==#
    # this gives distribution at the start for period t before period t shocks
    # have been realized
    omega1, omega1hat = Tfunc_general(omega_l, polp_l_check, pollamb_l, p, p.np, infl, Zl)
    pdist = sum(omega1, dims=2)

    # get implied aggregate Y
    # p_infl = p.pgrid .* exp(infl)
    Yimplied = sum((p.pgrid  .^ (-p.ϵ))' * pdist)

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0
    F = 0
    for pidx = 1:p.np
        for aidx = 1:p.na
            pval = p.pgrid[pidx]
            # pval = pval * exp(infl)
            aval = log(Z) .+ p.agrid[aidx]
            F += p.κ * pollamb[pidx, aidx] .* omega1hat[pidx, aidx] # who adjusts in a period
            Ld += pval^(-p.ϵ) * exp(-aval) * Yimplied * omega1[pidx,aidx]
        end
    end

    # monetary policy
    # talor rule
    r_val = rl +  p.ϕ_infl*(infl - p.Π_star) + p.ϕ_output*(Y-yss) + ϵ[2]
    mon_pol_error = r_val - r
    # Zmonerror = Zmon - p.ρ_agg * Zmonl - ϵ[2]

    cerror = C - Yimplied + F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = (1/Cl) - (1+r)*p.β*1/(exp(infl)*C) - η_ee
    zerror = log(Z) - p.ρ_agg * log(Zl) - ϵ[1]


    # == residual vector == #
    residvector = zero(X)
    residvector[1:sizedist] = vec(omega1 - omega)
    residvector[(sizedist+1):(2*sizedist)] = vec(V_l_check - V_l)


    # other error
    residvector[(2*sizedist+1):end] = [w-w_implied,Y-Yimplied,
                                        euler_error,cerror,
                                        mon_pol_error,
                                        zerror]

     
    return residvector
    
end

#==
Residual equations which equal 0 at equilibrium
Made for linear time iterations
Xl is lagged value of variables
X is current value of variables
Xf is future (expected) value of variables
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
function residequations_lti(Xl, X, Xf, ϵ, p)
    
    # default params
    @unpack np, na = p

    sizedist = np * na

    # unpacking
    omega_l = reshape(Xl[1:sizedist], np, na)
    omega = reshape(X[1:sizedist], np, na)
    omegaf = reshape(Xf[1:sizedist], np, na)

    V_l = reshape(Xl[(sizedist+1):(2*sizedist)], np, na)
    V = reshape(X[(sizedist+1):(2*sizedist)], np, na)
    V_f = reshape(Xf[(sizedist+1):(2*sizedist)], np, na)

    # need to rewrite this opening to accomodate vectors for V, polp, etc
    wl, rl, Yl, Cl, Zl  = exp.(Xl[(2*sizedist+1):(end-2)])
    w, r, Y, C, Z = exp.(X[(2*sizedist+1):(end-2)])
    wf, rf, Yf, Cf, Zf  = exp.(Xf[(2*sizedist+1):(end-2)])
    Zmonl, infl_l = Xl[(end-1):end]
    Zmon, infl = X[(end-1):end]
    Zmonf, infl_f = Xf[(end-1):end]

    stochdiscfactor = C/Cf

    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    # calculate implied polp
    V_check, Vadj, Vnoadj, polp_check, pollamb, _, _ = vBackwardFirm(
        (Y=Y, w=w), p, Z, V_f, infl, infl_f, stochdiscfactor = stochdiscfactor
    )

    #==
    Compute Distribution checks
    ==#
    omega1, omega1hat = Tfunc(omega_l, polp_check, pollamb, p)
    pdist = sum(omega1, dims=2)

    # get implied aggregate Y
    Yimplied = sum((p.pgrid .^ (-p.ϵ))' * pdist)

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0
    F = 0
    for pidx = 1:p.np
        for aidx = 1:p.na
            pval = p.pgrid[pidx]
            aval = Zl * p.agrid[aidx]
            F += p.κ * pollamb[pidx, aidx] .* omega1hat[pidx, aidx] # who adjusts in a period
            Ld += pval^(-p.ϵ) * exp(-aval) * Yl * omega1[pidx,aidx]
        end
    end

    # monetary policy
    mon_pol_error = infl - log(Yl) + log(Y) - Zmon
    Zmonerror = Zmon - p.ρ_agg * Zmonl - ϵ[2]

    cerror = C - Yimplied - F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = (1/C) - (1+rf)*p.β*1/(Cf)
    zerror = log(Z) - p.ρ_agg * log(Zl) - ϵ[1]

    # == residual vector == #
    residvector = zero(X)
    residvector[1:sizedist] = vec(omega1 - omega)
    residvector[(sizedist+1):(2*sizedist)] = vec(V_check - V)


    # other error
    residvector[(2*sizedist+1):end] = [w-w_implied,Y-Yimplied,
                                        euler_error,cerror,
                                        mon_pol_error,
                                        Zmonerror,
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