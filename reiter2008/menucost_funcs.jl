using QuantEcon, LinearAlgebra, Roots, Parameters
using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff

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
==#
function  vBackwardFirm(agg, params, Z::AbstractVector{T}, v1::AbstractVector{T}; 
                        stochdiscfactor = 1.0,
                        maxiter=10000, tol=1e-6, 
                        printinterval=1000, printinfo=true) where T


    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params
    

    # preallocate profit matrix
    profit_mat = eltype(Z)zeros(params.np, params.na);
    for pidx=1:np
        for aidx=1:na
            pval = params.pgrid[pidx];
            aval = Z * params.agrid[aidx];
            profit_mat[pidx, aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    error = 10;
    iter = 0;


    ev1 = v1 * aP';

    # iterate over choices
    vnoadjust = profit_mat + β * stochdiscfactor * ev1
    vadjust_val = profit_mat .- κ +  β * stochdiscfactor * ev1
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
    omega1hat = zeros(params.np, params.na);
    for pidx = 1:params.np
        for aidx = 1:params.na

            for a0idx = 1:params.na
                omega1hat[pidx, aidx] = omega1hat[pidx, aidx] + omega0[pidx, a0idx] * aP[a0idx, aidx];
            end
        end
    end

    # update policies
    omega1 = zeros(params.np, params.na);
    for pidx = 1:params.np
        for aidx = 1:params.na
            
            p1idx = polp[pidx, aidx];

            # non adjusters
            omega1[pidx, aidx] = omega1[pidx, aidx] + (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
            
            # adjusters
            omega1[p1idx, aidx] = omega1[p1idx, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];
        end
    end

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
function Tfunc_general(omega0, polp, pollamb, params, ngrid)
    
    aP = params.aP

    # update shock dist
    omega1hat = zeros(params.np, params.na);
    for pidx = 1:ngrid
        for aidx = 1:params.na

            for a0idx = 1:params.na
                omega1hat[pidx, aidx] = omega1hat[pidx, aidx] + omega0[pidx, a0idx] * aP[a0idx, aidx];
            end
        end
    end

    # update policies
    omega1 = zeros(ngrid, params.na);
    for pidx = 1:ngrid
        for aidx = 1:params.na
            
            pval = polp[aidx];

            # non adjusters
            omega1[pidx, aidx] = omega1[pidx, aidx] + (!pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            if pval > params.pgrid[1] && pval < params.pgrid[end]
                pidx_vals = searchsorted(params.pgrid, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid[pidx_hi] - params.pgrid[pidx_lo]

                wt_lo = (pval - params.pgrid[pidx_lo])/total_dist
                wt_hi = (params.pgrid[pidx_hi] - pval)/total_dist
                
                # adjusters
                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * pollamb[pidx, aidx] * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval <= params.pgrid[1]

                # adjusters
                omega1[1, aidx] = omega1[1, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval >= params.pgrid[end]

                # adjusters
                omega1[end, aidx] = omega1[end, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            end
        end
    end

    return omega1, omega1hat 
end

function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)


    omega1 = ones(params.np, params.na);
    omega1 = omega1 ./ (params.np*params.na);
    omega1hat = zeros(params.np, params.na);
    error = 10;
    iter = 0;
    
    while (error > tol) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc(omega0, polp, pollamb, params)
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
function findEquilibrium_ss(p; winit=1, tol=1e-4, max_iter=100, deltaw=0.1,
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
        Ld = 0
        F = 0
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
function residequations(Xl::AbstractArray{T}, X::AbstractVector{T}, 
                η::AbstractVector{T}, ϵ::AbstractVector{T}, 
                p)
    
    # default params
    @unpack np, na = p

    sizedist = np * na

    # unpacking
    omegahat_l = reshape(Xl[1:sizedist], np, na)
    omegahat = reshape(X[1:sizedist], np, na)
    omega_l = reshape(Xl[(sizedist+1):(2*sizedist)], np, na)
    omega = reshape(X[(sizedist+1):(2*sizedist)], np, na)

    Vadj_l = reshape(Xl[(2*sizedist+1):(3*sizedist)], np, na)
    Vadj = reshape(X[(2*sizedist+1):(3*sizedist)], np, na)

    Vnoadj_l = reshape(Xl[(3*sizedist+1):(4*sizedist)], np, na)
    Vnoadj = reshape(X[(3*sizedist+1):(4*sizedist)], np, na)
    @show eltype(Vadj_l)


    # need to rewrite this opening to accomodate vectors for V, polp, etc
    wl, rl, Yl, Cl, Zl  = exp.(Xl[(4*sizedist+1):(end-1)])
    w, r, Y, C, Z  = exp.(X[(4*sizedist+1):(end-1)])
    infl_l = Xl[end]
    infl = X[end]

    # 2 shocks - tfp and monetary
    ϵ_tfp = ϵ[1]
    ϵ_ngdp = ϵ[2]

    # expectation errors
    ηvadj = reshape(η[1:sizedist], np, na)
    ηvnoadj = reshape(η[(sizedist+1):(2*sizedist)], np, na)
    η_stoch, η_ee = η[(2*sizedist + 1):end]

    stochdiscfactor = Cl/C + η_stoch

    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    Vout = max.(Vadj, Vnoadj)
    # Vout_int = cubic_spline_interpolation((pgrid_orig, agrid_orig), Vout)
    # Vadj_l_check = zeros(p.np, p.na)
    # flowprofit = zeros(p.na)
    # for aidx=1:p.na
    #     pval = polp_l_val[aidx];
    #     aval = agrid_l[aidx];
    #     flowprofit[aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(wl/exp(aval))) * Yl;
    # end

    # calculate implied polp
    vout_l, Vadj_l_check, Vnadj_l_check, polp_l_check, pollamb_l, _, _ = vBackwardFirm(
        (Y=Yl, w=wl), p, Zl, Vout, stochdiscfactor = stochdiscfactor
    )
    polp_l_val_check = p.pgrid[polp_l_check]


    # store v adjust check
    # Vadj_l_check = zeros(p.np, p.na)
    # for pidx in 1:p.np
    #     for aidx in 1:p.na
    #         pval = polp_l_val[pidx]
    #         # compute expectations
    #         Ex = 0
    #         for a1idx in 1:p.na
    #             a1val = agrid[a1idx]
    #             Ex += Vout_int(pval, a1val) * p.aP[aidx, a1idx]
    #         end
    #         Vadj_l_check[pidx, aidx] = flowprofit[aidx] + β*stochdiscfactor*Ex
    #     end
    # end
    Vadj_l_check += ηvadj

    # store v noadjust check
    # Vnadj_l_check = zeros(p.np, p.na)
    # for pidx in 1:p.np
    #     for aidx in 1:p.na
    #         pval = p.pgrid[pidx]
    #         aval = p.agrid[aidx]
    #         val= (pval^(1-ϵ) - pval^(-ϵ)*(wl/exp(aval))) * Yl
    #         # compute expectations
    #         Ex = 0
    #         for a1idx in 1:p.na
    #             a1val = agrid[a1idx]
    #             Ex += Vout_int(pval, a1val) * p.aP[aidx, a1idx]
    #         end
    #         Vnadj_l_check[pidx, aidx] = val + β*stochdiscfactor*Ex
    #     end
    # end
    Vnadj_l_check += ηvnoadj

    pollamb = Vadj_l .> Vnoadj_l

    #==
    Compute Distribution checks
    ==#
    omega1, omega1hat = Tfunc(omega_l, polp_l_check, pollamb, p)
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
            F += p.κ * pollamb[pidx, aidx] .* omegahat[pidx, aidx] # who adjusts in a period
            Ld += pval^(-p.ϵ) * exp(-aval) * Yl * omega[pidx,aidx]
        end
    end

    # monetary policy
    mon_pol_error = infl - log(Yl) + log(Y) - ϵ[2]

    cerror = C - Yimplied - F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = (1/Cl) - (1+r)*p.β*1/(C) + η_ee
    zerror = log(Z) - p.ρ_agg * log(Zl) - ϵ[1]


    # == residual vector == #
    residvector = zero(X)
    residvector[1:sizedist] = vec(omega1hat - omegahat)
    residvector[(sizedist+1):(2*sizedist)] = vec(omega1 - omega)
    residvector[(2*sizedist+1):(3*sizedist)] = vec(Vadj_l_check - Vadj_l)
    residvector[(3*sizedist+1):(4*sizedist)] = vec(Vnadj_l_check - Vnoadj_l)


    # other error
    residvector[(4*sizedist+1):end] = [wl-w_implied,Yl-Yimplied,
                                        euler_error,cerror,
                                        mon_pol_error,
                                        zerror]
     
    return residvector
    
end

#==
Linarizeed coefficients
==#
function linearized_coeffs(equations, xss, shocks_sd, fargs)

    shocks_sd = atleast_2d(shocks_sd)
    shocks_ss = zero(shocks_sd)

    H1 = ForwardDiff.jacobian(t -> equations(t, xss,  shocks_ss, shocks_sd, fargs...), xss)
    H2 = ForwardDiff.jacobian(t -> equations(xss, t,  shocks_ss, shocks_sd, fargs...), xss)
    H3 = ForwardDiff.jacobian(t -> equations(xss, xss, t,  shocks_sd, fargs...), xss)
    H4 = ForwardDiff.jacobian(t -> equations(xss, xss, xss, t,  fargs...), shocks_ss)

    return H1, H2, H3, H4

end