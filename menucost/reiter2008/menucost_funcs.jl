using QuantEcon, LinearAlgebra, Roots, Parameters
# using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using SparseArrays, SparsityDetection
using Optim
using BSplineKit
using ApproxFun
using FastGaussQuadrature

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
Bellman operator for not adjusting
==#
function T_noadjust!(V, Vadj1, Vnoadj1, params ,agg; sdf=1.0, infl = 0.0)

    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params

    @inbounds for pidx = 1:np
        @inbounds for aidx = 1:na

            pval = pgrid[pidx] / (1.0 + infl)
            aval = agg.A .+ agrid[aidx]
            profit = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y

            Ex = 0.0
            @inbounds for a1idx = 1:na
                vnoadj1_fun = extrapolate(interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
                vnext = max(vnoadj1_fun(pval), Vadj1[a1idx])
                Ex += aP[aidx, a1idx] * vnext
            end
            
            valnoadjust = profit + β * sdf * Ex

            V[pidx, aidx] = valnoadjust

        end
    end


end

#==
Bellman operator for adjusting
    - choose optimal policy
    - overwrites polp and v
==#
function T_adjust!(polp, Vadj1, Vnoadj1, params ,agg; sdf=1.0)

    @unpack np, na, pgrid, agrid, ϵ, β, aP, σ, κ, plo, phi, Sad, Snoadj, xgh, wgh, ngh = params

    Vadj1_f = Fun(Sadj, ApproxFun.transform(Sadj, Vadj1));
    Vnoadj1_f = Fun(Snoadj, ApproxFun.transform(Snoadj, Vadj1));

    @inbounds for aidx = 1:na

        aval = agg.A .+ agrid[aidx]

        function objective(p1)
            profit = (p1^(1-ϵ) - p1^(-ϵ)*(agg.w/exp(aval))) * agg.Y - κ

            Ex = 0.0
            @inbounds for epsidx = 1:ngh
                a1val = aval + σ * xgh[epsidx]
                vnext = max(Vnoadj1_f(p1, a1val), Vadj1_f(a1val)) * wgh
                Ex += vnext
            end

            valadjust =  profit + β * sdf * Ex

            return -valadjust

        end
        
        result = optimize(objective, plo, phi, Brent())
        pstar = Optim.minimizer(result)
        polp[aidx] = pstar

        profit = (pstar^(1-ϵ) - pstar^(-ϵ)*(agg.w/exp(aval))) * agg.Y - κ

        Ex = 0.0
        @inbounds for a1idx = 1:na
            vnoadj1_fun = extrapolate(interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
            vnext = max(vnoadj1_fun(pstar), Vadj1[a1idx])
            Ex += aP[aidx, a1idx] * vnext
        end
        
        valadjust = profit + β * sdf * Ex

        V[aidx] = valadjust

    end


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
                    maxiter=1000, tol=1e-6, 
                    printinterval=50, printinfo=true,
                    howarditer=100)


    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ, plo, phi = params

    # initial values for v
    Vadjust0 = zeros(na)
    Vnoadjust0 = zeros(np, na)
    Vadjust1 = zeros(na)
    Vnoadjust1 = zeros(np, na)
    polp0 = collect(range(phi, plo, length=na))
    pollamb0 = BitArray(undef, np, na)

    error = 10;
    iter = 0;

    while (error > tol) && (iter < maxiter)

        for hidx=1:howarditer
            T_adjust_given!(Vadjust1, polp0, Vadjust0, Vnoadjust0, params, agg)
            T_noadjust!(Vnoadjust1, Vadjust0, Vnoadjust0, params, agg, infl=params.Π_star)
            Vadjust0 = deepcopy(Vadjust1)
            Vnoadjust0 = deepcopy(Vnoadjust1)
        end

        # iterate over choices
        polp1 = zero(polp0)
        T_adjust_max!(Vadjust1, polp1, Vadjust0, Vnoadjust0, params, agg)
        T_noadjust!(Vnoadjust1, Vadjust0, Vnoadjust0, params, agg, infl=params.Π_star)

        pollamb1 = repeat(Vadjust1', np, 1)  .> Vnoadjust1

        # error_adj = maximum(abs.(Vadjust1 - Vadjust0))
        # error_noadj = maximum(abs.(Vnoadjust1 - Vnoadjust0))
        # error = max(error_adj, error_noadj)
        error_polp = maximum(abs.(polp1 - polp0))
        error_pollamb = maximum(abs.(pollamb1 - pollamb0))
        error = max(error_polp, error_pollamb)

        polp0 = deepcopy(polp1)
        pollamb0 = deepcopy(pollamb1)

        iter += 1
        # Vadjust0 = deepcopy(Vadjust1)
        # Vnoadjust0 = deepcopy(Vnoadjust1)

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
   
    pollamb = repeat(Vadjust1', np, 1)  .> Vnoadjust1
    vout = max.(repeat(Vadjust1', np, 1), Vnoadjust1)

    return vout, Vadjust1, Vnoadjust1, polp0, pollamb0, iter, error

end


#==
Make policy functions denser in the p dimension
- only outputs dense pollab funcion
==#
function makedense(Vadjust, Vnoadjust, params, agg)

    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ, npdense, pgrid_dense = params

    pollamb_dense = zeros(npdense, na)

    for pidx = 1:npdense
        for aidx = 1:na
            pval = pgrid_dense[pidx]
            Vnoadjust_interp = interpolate(pgrid, Vnoadjust[:, aidx], BSplineOrder(4))
            pollamb_dense[pidx, aidx] = Vadjust[aidx] > Vnoadjust_interp(pval)
        end
    end

    return pollamb_dense
end

#==
yougn simulation of updating exogenous shock process
==#
function dist_updateshocks(omega0, params, Z)

    aP = params.aP

    omega1hat = zeros(params.npdense, params.na);
    agrid =  params.agrid
    for pidx = 1:params.npdense
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
function Tfunc_general(omega0, polp, pollamb, params, infl, Z)
    

    omega1hat = dist_updateshocks(omega0, params, Z)

    # update policies
    omega1 = zeros(params.npdense, params.na);
    

    for pidx = 1:params.npdense
        for aidx = 1:params.na

            # non adjusters
            pval0 = params.pgrid_dense[pidx]
            pval = pval0 / (1.0 + infl)
            aval =  params.agrid[aidx]
            if pval > params.plo && pval < params.phi
                pidx_vals = searchsorted(params.pgrid_dense, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid_dense[pidx_hi] - params.pgrid_dense[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid_dense[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1 - wt_lo

                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                

            elseif pval <= params.plo

                omega1[1, aidx] = omega1[1, aidx] +  (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            elseif pval >= params.phi

                omega1[end, aidx] = omega1[end, aidx] +  (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            end

            # adjusters
            pval0 = params.pgrid_dense[pidx]
            pval = polp[aidx]
            if pval > params.plo && pval < params.phi
                pidx_vals = searchsorted(params.pgrid_dense, pval)
                pidx_lo = last(pidx_vals)
                pidx_hi = pidx_lo + 1
                total_dist = params.pgrid_dense[pidx_hi] - params.pgrid_dense[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid_dense[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo
                
                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * pollamb[pidx, aidx] * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval <= params.plo

                omega1[1, aidx] = omega1[1, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pval >= params.phi

                omega1[end, aidx] = omega1[end, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            end
        end
    end

    # omega1 = omega1 ./ sum(omega1)
    # omega1hat = omega1hat ./ sum(omega1hat)
    return omega1, omega1hat

end


function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)


    omega1 = ones(params.npdense, params.na);
    omega1 = omega1 ./ (params.npdense*params.na);
    omega1hat = zeros(params.npdense, params.na);
    # omega1 = sparse(omega1)
    # omega1hat = sparse(omega1hat)
    error = 10;
    iter = 0;
    
    while (error > tol) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc_general(omega0, polp, pollamb, params, params.Π_star, 1.0)
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
x = [w , y]
==#
function equilibriumResidual(x, p)

    w = x[1]
    Y = x[2]
    agg = (w=w, Y=Y, A=0.0);
    v1, Vadjust, Vnoadjust, polp, pollamb, _, _  = viterFirm(agg, p; maxiter=10000, tol=1e-6, printinfo=false)
    pollamb_dense = makedense(Vadjust, Vnoadjust, params, agg)

    # get joint distribution of prices and shocks
    omegahat, omega = genJointDist(polp, pollamb_dense, p; printinfo=false);

    # get implied aggregate price
    aggprice = 0.0
    for pidx=1:p.npdense
        for aidx = 1:p.na
            pval = p.pgrid_dense[pidx]
            aggprice += pval^(1.0 - p.ϵ) * omega[pidx, aidx]
        end
    end
    aggprice = aggprice^(1.0/(1.0-p.ϵ))

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0.0
    F = 0.0
    for pidx = 1:p.npdense
        for aidx = 1:p.na
            pval = p.pgrid_dense[pidx]
            a1val = p.agrid[aidx]
            pchange = pollamb_dense[pidx, aidx]
            F += p.κ * pchange * omegahat[pidx, aidx]
            Ld += pval^(-p.ϵ) * exp(-a1val) * Y * omega[pidx,aidx]
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
    @unpack np, na, npdense, pgrid, pgrid_dense, aP, β, agrid = p

    sizedist = npdense * na
    sizev = np * na

    # unpacking
    omega_l = reshape(Xl[1:sizedist], npdense, na)
    omega = reshape(X[1:sizedist], npdense, na)

    Vadj_l = Xl[(sizedist+1):(sizedist+na)]
    Vadj = X[(sizedist+1):(sizedist+na)]

    Vnoadj_l = reshape(Xl[(sizedist+na+1):(sizedist+na+sizev) ], np, na)
    Vnoadj = reshape(X[(sizedist+na+1):(sizedist+na+sizev) ], np, na)

    polp_l = Xl[(sizedist+na+sizev+1):(sizedist+na+sizev+na)]
    polp = X[(sizedist+na+sizev+1):(sizedist+na+sizev+na)]

    wl, rl, Yl, Cl, Zl  = Xl[(sizedist+na+sizev+na+1):(end-1)]
    w, r, Y, C, Z = X[(sizedist+na+sizev+na+1):(end-1)]
    infl = Xl[end]
    # infl = X[end]

    # expectation errors
    ηv_adj = η[1:na]
    ηv_noadj = reshape(η[(na+1):(na+sizev)], np, na)
    η_polp = η[(na+sizev+1):(2*na+sizev)]
    # ηv_noadj = reshape(η[(na+1):(na + sizedist)], np, na)
    η_ee = η[2*na+sizev+1]


    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    sdf = Cl/C
    # calculate Vadjl, Vnoadjl given policy function
    Vadj_l_check = zero(Vadj)
    Vnoadj_l_check = zero(Vnoadj)
    V1 = max.(repeat(Vadj', np, 1), Vnoadj)
    agg = (Y=Y, w=w, A=Z)
    T_adjust_given!(Vadj_l_check, polp_l, Vadj, Vnoadj, params ,agg; sdf=sdf)
    T_noadjust!(Vnoadj_l_check, Vadj, Vnoadj, params ,agg; sdf=sdf, infl=infl)

    Vadj_l_check += ηv_adj
    Vnoadj_l_check += ηv_noadj

    # polp check
    polp_error = zeros(eltype(X), na)
    for aidx = 1:na
        aval = agrid[aidx]
        polpstar = polp_l[aidx]
        Ex = 0.0
        for a1idx = 1:na
            vprimeinter = Derivative(1) * extrapolate(interpolate(pgrid, V1[:, a1idx], BSplineOrder(4)), Smooth())
            Ex +=  aP[aidx, a1idx] * vprimeinter(polpstar)
        end
        dprofit = ((1.0 - p.ϵ)*polpstar^(-p.ϵ) + p.ϵ*polpstar^(-p.ϵ-1.0)*(w/exp(Z+aval))) * Y
        polp_error[aidx] = dprofit + β*sdf*Ex + η_polp[aidx]
    end

    #==
    Compute Distribution checks
    ==#
    # this gives distribution at the start for period t before period t shocks
    # have been realized
    pollamb_l_dense = makedense(Vadj_l, Vnoadj_l, p, agg)
    pollamb_dense = makedense(Vadj, Vnoadj, p, agg)

    omega1, omega0hat = Tfunc_general(omega_l, polp_l, pollamb_l_dense, p, infl, Z)
    # omega1hat = Tfunc_updateshocks(omega1, p, p.np_fine, Z)

    # get implied aggregate price
    aggprice = 0.0
    for pidx=1:p.npdense
        for aidx = 1:p.na
            pval = p.pgrid_dense[pidx]
            aggprice += pval^(1.0 - p.ϵ) * omega[pidx, aidx]
        end
    end
    aggprice = aggprice^(1.0/(1.0-p.ϵ))

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0.0
    F = 0.0
    for pidx = 1:p.npdense
        for aidx = 1:p.na
            pval = p.pgrid_dense[pidx]
            aval = p.agrid[aidx]
            pchange = pollamb_l_dense[pidx, aidx]

            F += p.κ * pchange * omega0hat[pidx, aidx]
            Ld += pval^(-p.ϵ) * exp(-aval) * Y * omega[pidx,aidx]
        end
    end

    # monetary policy
    # talor rule
    r_val = rl +  p.ϕ_infl*(infl - p.Π_star) + p.ϕ_output*(Y-yss) + ϵ[2]
    mon_pol_error = r_val - r

    cerror = C - Y + F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = 1.0/Cl - (1+rl)*p.β/(C*(1.0+infl)) - η_ee
    zerror = Z - p.ρ_agg * Zl - ϵ[1]


    # == residual vector == #
    residvector = zero(X)
    residvector[1:sizedist] = vec(omega1 - omega)
    residvector[(sizedist+1):(sizedist+na)] = vec(Vadj_l_check - Vadj_l)
    residvector[(sizedist+na+1):(sizedist+na+sizev) ] = vec(Vnoadj_l_check - Vnoadj_l)
    residvector[(sizedist+na+sizev+1):(sizedist+na+sizev+na)] = polp_error


    # other error
    residvector[(sizedist+na+sizev+na+1):end] = [w-w_implied,1.0 - aggprice,
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
