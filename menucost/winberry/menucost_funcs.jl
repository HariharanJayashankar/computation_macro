using QuantEcon, LinearAlgebra, Roots, Parameters
using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using Optim
using FastGaussQuadrature

#==
 This function maps z in [-1,1] into x in [xmin,xmax]; inverse of scale_down
 
 Inputs
   (1) z: array with elements in [-1,1]; if outside range, will force
   (2) xmin: scalar lower bound of range
   (3) xmax: scalar upper bound of range

 Outputs
   (1) x: array with elements in [xmin,xmax]
 
 Thomas Winberry, February 14th, 2018
==#
function scaleUp(z,xmin,xmax)

    # x = min(max((.5 .* (z .+ 1.0) .* (xmax - xmin)) .+ xmin,xmin .* ones(size(z))), 
    #     xmax .* ones(size(z)))

    x = (xmax - xmin) .* z .+ xmin .+ xmax
    x = x/2.0
    return x
end

#==
This function maps x in [xmin,xmax] into z in [-1,1]; inverse of scale_up

Inputs
  (1) x: array with elements in [xmin,xmax]; if outside range, will force
  (2) xmin: scalar lower bound of domain 
  (3) xmax: scalar upper bound of domain

Outputs
  (1) z: array with elements in [-1,1]

Thomas Winberry, Feburary 14th, 2018
==#
function scaleDown(x,xmin,xmax)

    # z = min(max(2.0 .* ((x .- xmin) / (xmax - xmin)) .- 1.0,-1.0 .* ones(size(x))),ones(size(x)));

    return z
end

#==
change of vars accroding to gauss legendre

if we want to calcuulate int_xlo^xhi x dx

We can do the transformation
x = ((xhi - xlo)*t + xhi + xlo)/2
=> t = (2x - xhi - xlo)/(xhi - xlo)
==#
function changeofvars(x, xlo, xhi)
    t = (2*x - xhi - xlo)/(xhi - xlo)
    return t
end

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
function  vBackwardFirm(agg, params, Z, Z_f, v1, infl, infl_f; 
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
            aval =  log(Z) + params.agrid[aidx];
            profit_mat[pidx, aidx] = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y;
        end
    end

    error = 10;
    iter = 0;

    # interpolate v1
    itp = interpolate(v1, (BSpline(Linear()), BSpline(Linear())))
    eitp = extrapolate(itp, Line())
    v1interp = Interpolations.scale(eitp, pgrid, agrid)
    v1_infl = zeros(eltype(v1), np, na)
    v1_noadj = zeros(eltype(v1), np, na)
    for pidx = 1:np
        for aidx = 1:na
            pval = pgrid[pidx]
            pval_adj = pval / exp(infl_f)

            aval = log(Z_f) + agrid[aidx]
            v1_infl[pidx, aidx] = v1interp(pval_adj, aval)
        end
    end

    ev1_infl = v1_infl * aP';

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

    return omega1, omega1hat

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
#==
For a given set of moments and distributino parameters
get the density for (a, p)
m and g are matrices of size (ng, ng), where
ng is the level of approximation

Notation:
m[i,j] is m_i^j in winberry notation
==#
function getDensity(p, a, m, g, g0, params)

    @unpack ng = params

    # first moments
    gout = g[1] * (p - m[1]) + g[2] * (a - m[2])

    # rest of the moments
    idx=2
    for i=2:ng
        for j=0:i
            idx += 1
            inner = (p - m[1])^(i-j) * (a - m[2])^j - m[idx]
            gout += g[idx] * inner
        end
    end

    gout = g0 * exp(gout)

    return gout

end


#==
Objective function for getting distributional parameters
given a set of moments
m is a matrix of moments of size (ng, ng)
g is passed in as a long vector which is then reshaped into a matrix
so it can be passed into getDensity

See written notes to see how vector is filled up into a matrix form
For example if ng = 2
We have the parameters {g_1^1, g_1^2, g_2^0, g_2^1, g_2^2}
The amtrix will be
[ x g_1^1 g_1^2
 g_2^0, g_2^2, g_2^2]
==#
function objectiveDensity(g, m, params)

    @unpack plo, phi, pgrid, agrid, np, na, ng, nquad, xquad, wquad = params

    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    pscale = (phi - plo)/2.0
    agrid_quad = scaleUp(xquad, alo, ahi)
    ascale = (ahi - alo)/2.0

    integral = 0.0
    for pidx = 1:nquad
        for aidx = 1:nquad
            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]
            integral += getDensity(pval, aval, m, g, 1.0, params) * wquad[pidx] * wquad[aidx]
        end
    end

    integral *= pscale * ascale

    return integral

end

#==
Get density gradient to help with optimizer
==#
function getDensityG!(G, x, m, params)


    @unpack plo, phi, pgrid, agrid, np, na, ng, nquad, xquad, wquad = params

    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    agrid_quad = scaleUp(xquad, alo, ahi)
    pscale = (phi - plo)/2.0
    ascale = (ahi - alo)/2.0

    # set to 0 initially
    G[:] .= 0.0

    # first two moments
    # need to integrate
    for pidx = 1:nquad
        for aidx = 1:nquad
            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]

            moment = pval - m[1]
            G[1] += moment * getDensity(pval, aval, m, x, 1.0, params) * wquad[pidx] * wquad[aidx]
        end
    end
    for pidx = 1:nquad
        for aidx = 1:nquad
            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]

            moment = aval - m[2]
            G[2] += moment * getDensity(pval, aval, m, x, 1.0, params) * wquad[pidx] * wquad[aidx]
        end
    end

    # rest of the moments
    idx=2
    for i=2:ng
        for j=0:i
            # need to integrate
            idx += 1
            for pidx = 1:nquad
                for aidx = 1:nquad

                    pval = pgrid_quad[pidx]
                    aval = agrid_quad[aidx]

                    moment = (pval - m[1])^(i-j) * (aval - m[2])^j - m[idx]
                    # need to offset G by 2 more since we dont divide into G_Base and G_rest
                    G[idx] += moment * getDensity(pval, aval, m, x, 1.0, params) * wquad[pidx] * wquad[aidx]
                end
            end
        end
    end

    G[:] *= pscale * ascale

end

#==
Iterate parametrized distribution by one period using policy rules,
the distribution parameters and the starting moments
==#
function iterateDist(g0, g, m0, polp, pollamb, params, infl, Zl, Z)

    @unpack np, na, plo, phi, pgrid, agrid, ng, 
        nquad, xquad, wquad, ngh, xgh, wgh,
        ρ, σ  = params
    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    pscale = (phi - plo)/2.0
    agrid_quad = scaleUp(xquad, alo, ahi)
    ascale = (phi - plo)/2.0

    m1 = zero(m0)

    # interpolate pol
    itp = interpolate(pgrid[polp], (BSpline(Linear()), BSpline(Linear())))
    eitp = extrapolate(itp, Line())
    polp_interp = Interpolations.scale(eitp, pgrid, agrid)

    itp = interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
    eitp = extrapolate(itp, Line())
    pollamb_interp = Interpolations.scale(eitp, pgrid, agrid)

    for pidx = 1:nquad
        for aidx = 1:nquad

            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]
            aval = log(Zl) + aval

            # iterate over next a
            for epsidx = 1:ngh
                epsval = xgh[epsidx]
                a1val = ρ*aval + σ*epsval
                a1val = log(Z) + a1val

                # remember a1s realize before next periods p
                p1val = polp_interp(pval, a1val)
                pchange = pollamb_interp(pval, a1val)

                p1val = pchange * p1val + (1.0 - pchange) * pval
                p1val = p1val / exp(infl)

                density = getDensity(pval, aval, m0, g, g0, params) * wquad[pidx] * wquad[aidx] * wgh[epsidx]

                m1[1] += density * p1val
                m1[2] += density * a1val

                # higher moments
                idx = 2
                for i=2:ng
                    for j=0:i
                        idx += 1
                        m1[idx] = ((p1val - m1[1])^(i-j)) * ((a1val - m1[2])^j) * density
                    end
                end


            end
        end
    end

    m1 *= pscale * ascale
    return m1

end


function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)

    @unpack ng, pgrid, agrid, dampening = params

    tol_hist = tol * 1e3
    # == initial histogram approach to get decent starting moments == #
    omega1 = ones(params.np, params.na);
    omega1 = omega1 ./ (params.np*params.na);
    omega1hat = zeros(params.np, params.na);
    error = 10;
    iter = 0;
    
    while (error > tol_hist) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc_general(omega0, polp, pollamb, params, params.np, params.Π_star, 1.0)
        error = maximum(abs.(omega1 - omega0))
        iter += 1;
        omega0hat = omega1hat;

        if ((iter==1) || (mod(iter, printinterval) == 0)) & printinfo
            println("Histogram iterations: $iter, Error $error")
        end

    end

    if printinfo
        println("Final Histogram Iterations: $iter, Final error: $error")
    end

    # == Winberry == #
    # calculate moments from histogram
    m0 = zeros(p.nparams)
    pdist = sum(omega1, dims=2)
    adist = sum(omega1, dims=1)
    m0[1] = sum(pdist .* pgrid)
    m0[2] = sum(adist .* agrid)


    idx = 2

    for i=2:ng
        for j=0:i
            inner = 0.0

            for pidx = 1:params.np
                for aidx = 1:params.na
                    pval = pgrid[pidx]
                    aval = agrid[aidx]
                    inner += (pval - m0[1])^(i-j) * (aval - m0[2])^j * omega1[pidx, aidx]
                end
            end

            idx += 1
            m0[idx] = inner
        end
    end

    # == iterate to steady state == #
    error = 10;
    iter = 0;

    # some initial guesses for g
    gprev = 0.2 * ones(p.nparams)

    while (error > tol) && (iter < maxiter)

        # find parameters
        result = optimize(
            x -> objectiveDensity(x, m0, p),
            (G,x) -> getDensityG!(G, x, m0, p),
            gprev,
            AdaMax(),
            Optim.Options(x_tol=1e-3, f_tol=1e-3, g_tol=1e-3, iterations=1_000_000)
        )
        gest = Optim.minimizer(result)
        densityOut = Optim.minimum(result)
        if Optim.converged(result) == false
            @error "Optimizer failed to get distributional parameters"
            @show m0
        end
        g0 = 1.0 / densityOut

        # iterate LOM
        m1 = iterateDist(g0, gest, m0, polp, pollamb, params, 0.0, 1.0, 1.0)
        error = maximum(abs.(m1 - m0))
        m0 = dampening * m1 + (1.0 - dampening) * m0
        iter += 1
        gprev = gest

        if ((iter==1) || (mod(iter, printinterval) == 0)) & printinfo
            println("Parametrized Dist iterations: $iter, Error $error")
        end


    end

    if printinfo
        println("Final Parametrized Dist Iterations: $iter, Final error: $error")
    end

    # get the parameters for final output
    result = optimize(
        x -> objectiveDensity(x, m0, p),
        (G,x) -> getDensityG!(G, x, m0, p),
        gprev,
        AdaMax(),
        Optim.Options(x_tol=1e-3, f_tol=1e-3, g_tol=1e-3, iterations=1_000_000)
    )
    gest = Optim.minimizer(result)
    densityOut = Optim.minimum(result)
    g0 = 1.0 / densityOut

    return m0, g0, gest


end

#==
Find equilibiurm Y and w to clear steady state
==#
function findEquilibrium_ss(p; winit=1, tol=1e-3, max_iter=200, deltaw=0.05,
                            Yinit=1, deltaY=0.05,
                            printinterval=10)
    
    @unpack np, na, plo, phi, pgrid, agrid, ng, 
        nquad, xquad, wquad, ngh, xgh, wgh, ρ, σ, nparams = p
    
    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    agrid_quad = scaleUp(xquad, alo, ahi)
    pscale = (phi - plo)/2.0
    ascale = (ahi - alo)/2.0


    # gettngg value funcion
    w0 = winit;
    Y0 = Yinit
    iter = 0;
    error = 10

    # preallocating
    Vadjust = zeros(p.np, p.na)
    Vnoadjust = zeros(p.np, p.na)
    polp = zeros(p.np, p.na)
    pollamb = zeros(p.np, p.na)
    moments = zeros(nparams)
    g = zeros(nparams)
    g0 = 0.0
    C = 0.0

    # outer labour market loop
    while iter < max_iter && error > tol

        agg = (w=w0, Y=Y0);
        V, Vadjust ,Vnoadjust, polp, pollamb  = viterFirm(agg, p; maxiter=10000, tol=1e-6, printinfo=false)

        # get joint distribution of prices and shocks
        moments, g0, g = genJointDist(polp, pollamb, p; printinfo=false);

        # get implied aggregate Y
        Yimplied = 0.0
        for pidx in 1:nquad
            for aidx in 1:nquad
                pval = pgrid_quad[pidx]
                aval = agrid_quad[aidx]
                density = getDensity(pval, aval, moments, g, g0, p) * wquad[pidx] * wquad[aidx]
                Yimplied += pval ^ (-p.ϵ) * density
            end
        end
        Yimplied *= pscale * ascale

        # get profits to give HH
        # get aggregate fixed cost payments
        # get labour demand
        # integrate
        # interpolate pol
        itp = interpolate(pgrid[polp], (BSpline(Linear()), BSpline(Linear())))
        eitp = extrapolate(itp, Line())
        polp_interp = Interpolations.scale(eitp, pgrid, agrid)

        itp = interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
        eitp = extrapolate(itp, Line())
        pollamb_interp = Interpolations.scale(eitp, pgrid, agrid)
        Ld = 0.0
        F = 0.0
        for pidx = 1:nquad
            for aidx = 1:nquad
                pval = pgrid_quad[pidx]
                aval = agrid_quad[aidx]
                density = getDensity(pval, aval, moments, g, g0, p) * wquad[pidx] * wquad[aidx]
                Ld += pval^(-p.ϵ) * exp(-aval) * Y0 * density

                # F is done after realizing shocks
                for epsidx = 1:ngh
                    epsval = xgh[epsidx]
                    a1val = ρ*aval + σ*epsval
                    F += p.κ * pollamb_interp(pval, a1val)*density*wgh[epsidx] # who adjusts in a period
                end
            end
        end
        F *= pscale * ascale
        Ld *= pscale * ascale

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

    return w0, Y0, Vadjust, Vnoadjust, polp, pollamb,  moments, g, g0, C, iter, error
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
    @unpack np, na, ng, nparams, agrid, pgrid, plo, phi, 
        xquad, wquad, nquad, ngh, xgh, wgh, ρ, σ = p

    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    pscale = (phi - plo)/2.0
    agrid_quad = scaleUp(xquad, alo, ahi)
    ascale = (ahi - alo)/2.0

    sizeval = np * na

    # unpacking distributions
    g0_l = Xl[1]
    g0 = X[1]
    g_l = Xl[2:(1+nparams)]
    g = X[2:(1+nparams)]

    m0 = Xl[(2+nparams):(2+2*nparams-1)]
    m1 = X[(2+nparams):(2+2*nparams-1)]

    # unpacking value functions
    Vadj_l = reshape(Xl[(2+2*nparams):(2+2*nparams+sizeval-1)], np, na)
    Vadj = reshape(X[(2+2*nparams):(2+2*nparams+sizeval-1)], np, na)
    Vnoadj_l = reshape(Xl[(2+2*nparams+sizeval):(1+2*nparams+2*sizeval)], np, na)
    Vnoadj = reshape(X[(2+2*nparams+sizeval):(1+2*nparams+2*sizeval)], np, na)

    polp_l = reshape(Xl[(2+2*nparams+sizeval):(1+2*nparams+2*sizeval)], np, na)
    polp = reshape(X[(2+2*nparams+sizeval):(1+2*nparams+2*sizeval)], np, na)

    wl, rl, Yl, Cl, Zl  = Xl[(2+2*nparams+2*sizeval):(end-1)]
    w, r, Y, C, Z = X[(2+2*nparams+2*sizeval):(end-1)]
    infl_l = Xl[end]
    infl = X[end]

    # expectation errors
    ηv = reshape(η[1:sizeval], np, na)
    η_ee = η[end]

    stochdiscfactor = (Cl*exp(infl))/C 

    #==
    Compute Value functions given optimal adjusting price rule
    ==#
    # calculate implied polp
    V = max.(Vadj, Vnoadj)
    V_l_check, Vadj_l_check, Vnoadj_l_check, polp_l_check, pollamb_l, _, _ = vBackwardFirm(
        (Y=Yl, w=wl), p, Zl, Z, V, infl_l, infl, stochdiscfactor = stochdiscfactor
    )
    V_l_check += ηv

    pollamb = Vadj .> Vnoadj

    #==
    Compute Distribution checks
    ==#
    # this gives distribution at the start for period t before period t shocks
    # have been realized
    m1_check = iterateDist(g0_l, g_l, m0, polp_l_check, pollamb_l, p, infl, Zl, Z)
    # get the parameters for final output
    result = optimize(
        x -> objectiveDensity(x, m1, p),
        (G,x) -> getDensityG!(G, x, m1, p),
        g_l,
        AdaMax(),
        Optim.Options(x_tol=1e-3, f_tol=1e-3, g_tol=1e-3, iterations=1_000_000)
    )
    
    gcheck = Optim.minimizer(result)
    densityOut = Optim.minimum(result)
    g0_check = 1.0 / densityOut

    #== 
    compute implied values from distribution
    ==#

    # get implied aggregate Y
    Yimplied = 0.0
    for pidx in 1:nquad
        for aidx in 1:nquad
            pval = pgrid_quad[pidx]
            aval = log(Z) + agrid_quad[aidx]
            density = getDensity(pval, aval, m1, g, g0, p) * wquad[pidx] * wquad[aidx]
            Yimplied += pval ^ (-p.ϵ) * density
        end
    end
    Yimplied *= pscale * ascale * Y

    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    # interpolate pol
    itp = interpolate(pgrid[polp], (BSpline(Linear()), BSpline(Linear())))
    eitp = extrapolate(itp, Line())
    polp_interp = Interpolations.scale(eitp, pgrid, agrid)

    itp = interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
    eitp = extrapolate(itp, Line())
    pollamb_interp = Interpolations.scale(eitp, pgrid, agrid)
    Ld = 0.0
    F = 0.0
    for pidx = 1:nquad
        for aidx = 1:nquad
            pval = pgrid_quad[pidx]
            aval = log(Zl) + agrid_quad[aidx]
            density = getDensity(pval, aval, m1, g, g0, p) * wquad[pidx] * wquad[aidx]

            # F is done after realizing shocks
            for epsidx = 1:ngh
                epsval = xgh[epsidx]
                a1val = log(Z) + ρ*aval + σ*epsval
                pchange =  pollamb_interp(pval, a1val)
                p1val = pchange * polp_interp(pval, a1val) + (1.0 - pchange)*pval
                F += p.κ * pchange*density*wgh[epsidx] # who adjusts in a period
                Ld += p1val^(-p.ϵ) * exp(-a1val) * Yl * density*wgh[epsidx]
            end
        end
    end
    F *= pscale * ascale
    Ld *= pscale * ascale

    # monetary policy
    # talor rule
    r_val = rl +  p.ϕ_infl*(infl - p.Π_star) + p.ϕ_output*(Y-yss) + ϵ[2]
    mon_pol_error = r_val - r

    cerror = C - Y + F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    euler_error  = (1/Cl) - (1+rl)*p.β*1/(exp(infl)*C) - η_ee
    zerror = log(Z) - p.ρ_agg * log(Zl) - ϵ[1]


    # == residual vector == #
    residvector = zero(X)

    # distributions
    residvector[1] = g0_check - g0
    residvector[2:(1+nparams)] = gcheck - g
    residvector[(2+nparams):(2+2*nparams-1)] = vec(m1_check) - vec(m1)

    # value function
    residvector[(2+2*nparams):(2+2*nparams+sizeval-1)] = vec(Vadj_l) - vec(Vadj_l_check)
    residvector[(2+2*nparams+sizeval):(1+2*nparams+2*sizeval)] = vec(Vnoadj_l) - vec(Vnoadj_l_check)


    # other errors
    residvector[(2+2*nparams+2*sizeval):end] = [w-w_implied,Y-Yimplied,
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
