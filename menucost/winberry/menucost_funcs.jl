using QuantEcon, LinearAlgebra, Roots, Parameters
using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using Optim
using FastGaussQuadrature
using BSplineKit

function findnearest(a,x)
       length(a) > 0 || return 0:-1
       r = searchsorted(a,x)
       length(r) > 0 && return r
       last(r) < 1 && return searchsorted(a,a[first(r)])
       first(r) > length(a) && return searchsorted(a,a[last(r)])
       x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
       x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
       return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end

@views function simps(f::Function, x::AbstractVector)
    n = length(x) - 1
    h = (x[end]-x[begin])/n
    I= h/3*(f(x[1])+2*sum(f,x[3:2:end-2])+4*sum(f,x[2:2:end])+f(x[end]))
    return I
end
simps(f::Function, a::Real, b::Real, n::Integer) = simps(f, a:((b-a)/n):b)
simps(x::AbstractVector) = simps(identity, x)

function simps2d(func::Function, a::Real, b::Real, c::Real, d::Real, NX::Integer, NY::Integer)
    
    # Ensure the number of intervals is even!
    NX = Int(2*ceil(NX/2));
    NY = Int(2*ceil(NY/2));
    # Set up the integration step sizes in the x and y directions
    hx = (b - a)/NX;
    hy = (d - c)/NY;
    # define grid vectors
    xg = a:hx:b;
    yg = c:hy:d;
    # xxg, yyg = meshgrid(xg,yg);
    # Now set up a matrix U that contains the values of the function evaluated at all
    # points on the 2D grid setup by xg and yg.
    U = func.(xg', yg);
    # Evaluate the contribution from the corner points first.
    # These all have weight 1. NB U(1,1) corresponds to func(a,b) etc.
    s1 = ( U[1,1] + U[1,NY+1] + U[NX+1,1] + U[NX+1,NY+1] );
    # Now sum the contributions from the terms along each edge not including
    # corners. There are 4 edges in the 2D case that contribute to the sum
    # and we have points with weight 4 and points with weight 2. Points
    # with weight 4 are acessed by indices 2:2:N (N=NX,NY,NZ), while points with
    # weight 2 are accessed by indices 3:2:N-1.
    # Define vectors of odd and even indices for each direction:
    ixo = 2:2:NX;
    ixe = 3:2:NX-1;
    iyo = 2:2:NY;
    iye = 3:2:NY-1;
    s2 = 2*( sum(U[1,iye]) + sum(U[NX+1,iye]) + sum(U[ixe,1]) + sum(U[ixe,NY+1]) );
    s3 = 4*( sum(U[1,iyo]) + sum(U[NX+1,iyo]) + sum(U[ixo,1]) + sum(U[ixo,NY+1]) );
    # Now we look at the remaining contributions on the interior grid points. 
    # Looking at our array example above we see that there
    # are only 3 different weights viz. 16, 8 and 4. Some thought will show that
    # using our definitions above for odd and even gridpoints, that weight 16 is
    # only found at points (xodd, yodd), weight 4 is found at points (xeven,yeven)
    # while weight 8 is found at both (xodd,yeven) or (xeven,yodd).
    # Our contribution from interior points is then
    s4 = 16*sum( sum( U[ixo,iyo] ) ) + 4*sum( sum( U[ixe,iye] ) );
    s5 =  8*sum( sum( U[ixe,iyo] ) ) + 8*sum( sum( U[ixo,iye] ) );
    # Finally add all the contributions and multiply by the step sizes hx, hy and 
    # a factor 1/9 (1/3 in each direction).
    out = s1 + s2 + s3 + s4 + s5;
    out = out*hx*hy/9.0

    return out

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
T_adjusting_given
- update value function for adjusting taking as given the policy
- overwrites V
==#
function T_adjust_given!(V, polp, Vadj1, Vnoadj1, params ,agg; sdf=1.0)

    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params

    @inbounds for aidx = 1:na
        pval = polp[aidx]
        aval = agg.A .+ agrid[aidx]

        profit = (pval^(1-ϵ) - pval^(-ϵ)*(agg.w/exp(aval))) * agg.Y - κ

        Ex = 0.0
        @inbounds for a1idx = 1:na
            vnoadj1_fun = BSplineKit.extrapolate(BSplineKit.interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
            vnext = max(vnoadj1_fun(pval), Vadj1[a1idx])
            Ex += aP[aidx, a1idx] * vnext
        end
        
        valadjust = profit + β * sdf * Ex

        V[aidx] = valadjust
        

    end


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
                vnoadj1_fun = BSplineKit.extrapolate(BSplineKit.interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
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
function T_adjust_max!(V, polp, Vadj1, Vnoadj1, params ,agg; sdf=1.0)

    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ, plo, phi = params

    @inbounds for aidx = 1:na

        aval = agg.A .+ agrid[aidx]

        function objective(p1)
            profit = (p1^(1-ϵ) - p1^(-ϵ)*(agg.w/exp(aval))) * agg.Y - κ

            Ex = 0.0
            @inbounds for a1idx = 1:na
                vnoadj1_fun = BSplineKit.extrapolate(BSplineKit.interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
                vnext = max(vnoadj1_fun(p1), Vadj1[a1idx])
                Ex += aP[aidx, a1idx] * vnext
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
            vnoadj1_fun = BSplineKit.extrapolate(BSplineKit.interpolate(pgrid, Vnoadj1[:, a1idx], BSplineOrder(4)), Smooth())
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
            Vnoadjust_interp = BSplineKit.interpolate(pgrid, Vnoadjust[:, aidx], BSplineOrder(4))
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
            pidx_vals = searchsorted(params.pgrid_dense, pval)
            pidx_lo = last(pidx_vals)
            pidx_hi = pidx_lo + 1
            if pidx_hi < params.npdense && pidx_lo > 1

                total_dist = params.pgrid_dense[pidx_hi] - params.pgrid_dense[pidx_lo]
                wt_lo = 1.0 - (pval - params.pgrid_dense[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1 - wt_lo

                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];
                

            elseif pidx_lo <= 1

                omega1[1, aidx] = omega1[1, aidx] +  (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            elseif pidx_hi >= params.npdense

                omega1[end, aidx] = omega1[end, aidx] +  (1.0 - pollamb[pidx, aidx]) * omega1hat[pidx, aidx];

            end

            # adjusters
            pval0 = params.pgrid_dense[pidx]
            pval = polp[aidx]
            pidx_vals = searchsorted(params.pgrid_dense, pval)
            pidx_lo = last(pidx_vals)
            pidx_hi = pidx_lo + 1
            if pidx_hi < params.npdense && pidx_lo > 1
                total_dist = params.pgrid_dense[pidx_hi] - params.pgrid_dense[pidx_lo]

                wt_lo = 1.0 - (pval - params.pgrid_dense[pidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo
                
                omega1[pidx_hi, aidx] = omega1[pidx_hi, aidx] + wt_hi * pollamb[pidx, aidx] * omega1hat[pidx, aidx];
                omega1[pidx_lo, aidx] = omega1[pidx_lo, aidx] + wt_lo * pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pidx_lo <= 1

                omega1[1, aidx] = omega1[1, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            elseif pidx_hi >= params.npdense

                omega1[end, aidx] = omega1[end, aidx] + pollamb[pidx, aidx] * omega1hat[pidx, aidx];

            end
        end
    end

    # omega1 = omega1 ./ sum(omega1)
    # omega1hat = omega1hat ./ sum(omega1hat)
    return omega1, omega1hat

end




#==
For a given set of moments and distributino parameters
get the density for (a, p)
m and g are matrices of size (ng, ng), where
ng is the level of approximation
==#
function getDensity(p, m, g, g0, params)

    @unpack ng = params

    # first moments
    gout = g[1] * (p - m[1])

    # rest of the moments
    for i=2:ng
        inner = (p - m[1])^(i) - m[i]
        gout += g[i] * inner
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

    @unpack plo, phi, pgrid, agrid, np, na, ng, nsimp = params

    # simpsons quadrature
    integral = simps(p -> getDensity(p, m, g, 1.0, params), plo, phi, nsimp)

    return integral

end

#==
Get density gradient to help with optimizer
==#
function getDensityG!(G, x, m, params)


    @unpack plo, phi, pgrid, agrid, np, na, ng, nsimp = params


    
    # simpsons quadrature
    # first two moments
    # need to integrate
    
    G[1] = simps(p -> (p - m[1]) * getDensity(p, m, x, 1.0, params), plo, phi, nsimp)

    # rest of the moments
    for i=2:ng
        # need to integrate
        G[i] = simps(p -> ((p - m[1])^i - m[i]) * getDensity(p, m, x, 1.0, params), plo, phi, nsimp)
    end

end



#==
Iterate parametrized distribution by one period using policy rules,
the distribution parameters and the starting moments
==#
function iterateDist(g0, g, m0, polp, pollamb, params, infl)

    @unpack np, na, plo, phi, pgrid, pgrid_dense, agrid, ng, 
        aP, ρ, σ, nsimp, aPstationary  = params

    m1 = zero(m0)

    polp_f = BSplineKit.interpolate(agrid, polp, BSplineOrder(4))
    polp_f = BSplineKit.extrapolate(polp_f, Smooth())


    # for each moment
    # first moment
    for a1idx = 1:na
        function inner_integral(pval, a1idx)
            mout = 0.0
            for aidx = 1:na

                aval = agrid[aidx]
                densityP = getDensity(pval, m0[:, aidx], g[:, aidx], g0[aidx], params)

                pval_nearest = findnearest(pgrid_dense, pval)[1]
                pchange = pollamb[pval_nearest, aidx]
                density = densityP * aP[aidx, a1idx]

                p1val = pchange * polp_f(aval) + (1.0 - pchange) * pval/(1.0 + infl)

                mout += p1val * density
            end
            return mout
        end

        m1[1, a1idx] = simps(p -> inner_integral(p, a1idx), plo, phi, nsimp)

        # for higher moments
        for j = 2:ng
            function inner_integral(pval,j,a1idx)
                mout = 0.0
                for aidx = 1:na

                    densityP = getDensity(pval, m0[:, aidx], g[:, aidx], g0[aidx], params)
                    aval = agrid[aidx]
                    pval_nearest = findnearest(pgrid_dense, pval)[1]
                    pchange = pollamb[pval_nearest, aidx]

                    density = densityP * aP[aidx, a1idx]

                    p1val = pchange * polp_f(aval) + (1.0 - pchange) * pval/(1.0 + infl)

                    mout += ((p1val - m1[1,a1idx])^j) * density

                end
                return mout
            end
            m1[j, a1idx] = simps(p -> inner_integral(p,j,a1idx), plo, phi, nsimp)

        end
    end

    return m1
end


function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)

    @unpack na, ng, pgrid, pgrid_dense, agrid, dampening = params

    tol_hist = tol
    # == initial histogram approach to get decent starting moments == #
    omega1 = ones(params.npdense, params.na);
    omega1 = omega1 ./ (params.npdense*params.na);
    omega1hat = zeros(params.npdense, params.na);
    error = 10;
    iter = 0;
    
    while (error > tol_hist) && (iter < maxiter)

        omega0 = omega1
        omega1, omega1hat = Tfunc_general(omega0, polp, pollamb, params, params.Π_star, 1.0)
        error = maximum(abs.(omega1 - omega0))
        iter += 1;
        omega0 = omega1;

        if ((iter==1) || (mod(iter, printinterval) == 0)) & printinfo
            println("Histogram iterations: $iter, Error $error")
        end

    end

    if printinfo
        println("Final Histogram Iterations: $iter, Final error: $error")
    end

    # == Winberry == #
    # calculate moments from histogram
    # m0 is a set of conditional moments conditional on each aidx
    m0 = zeros(ng, na)
    for aidx = 1:na
        pdist = omega1[:, aidx]
        pdist = pdist ./ sum(pdist) # sum to one
        m0[1, aidx] = sum(pdist .* pgrid_dense)

        for i=2:ng
            inner = 0.0

            for pidx = 1:params.npdense
                pval = pgrid_dense[pidx]
                inner += (pval - m0[1, aidx])^(i) * pdist[pidx]
            end

            m0[i, aidx] = inner
        end

    end

    # == iterate to steady state == #
    error = 10;
    iter = 0;
    # m0 = ones(ng, na) * p.pflex
    @show m0

    # some initial guesses for g
    # gprev = zeros(ng, na)
    gprev=  ones(ng, na) * -1e-2
    g0 = zeros(na)
    while (error > tol) && (iter < maxiter)

        # find parameters
        m1 = zero(m0)
        g1 = zero(m0)
        for aidx = 1:na
            result = optimize(
                x -> objectiveDensity(x, m0[:, aidx], p),
                (G,x) -> getDensityG!(G, x, m0[:, aidx], p),
                gprev[:, aidx],
                LBFGS()
            )
            gprev[:,aidx] = Optim.minimizer(result)
            densityOut = Optim.minimum(result)
            if Optim.converged(result) == false
                @show m0
                @error "Optimizer failed to get distributional parameters"
            end
            g0[aidx] = 1.0 / densityOut
            @show m0
            @show g0
        end

        m1 = iterateDist(g0, gprev, m0, polp, pollamb, params, 0.0)
        g1 = gprev

        error = maximum(abs.(m1 - m0))
        m0 = dampening * m1 + (1.0 - dampening) * m0
        iter += 1
        gprev = g1

        if ((iter==1) || (mod(iter, printinterval) == 0)) & printinfo
            println("Parametrized Dist iterations: $iter, Error $error")
        end


    end

    if printinfo
        println("Final Parametrized Dist Iterations: $iter, Final error: $error")
    end

    # get the parameters for final output
    gout = zero(m0)
    g0out = zeros(na)
    for aidx = 1:na
        result = optimize(
            x -> objectiveDensity(x, m0[:, aidx], p),
            (G,x) -> getDensityG!(G, x, m0[:, aidx], p),
            gprev[:, aidx],
            LBFGS()
        )
        gest = Optim.minimizer(result)
        densityOut = Optim.minimum(result)
        g0 = 1.0 / densityOut

        gout[:, aidx] = gest
        g0out[aidx] = g0
    end

    return m0, g0out, gout, omega1, omega1hat


end

#==
x = [w , y]
==#
function equilibriumResidual(x, p)

    @unpack np, na, plo, phi, pgrid, agrid, ng, 
        ρ, σ = p
    
    alo = minimum(agrid)
    ahi = maximum(agrid)
    pgrid_quad = scaleUp(xquad, plo, phi)
    agrid_quad = scaleUp(xquad, alo, ahi)
    pscale = (phi - plo)/2.0
    ascale = (ahi - alo)/2.0

    w = x[1]
    Y = x[2]

    agg = (w=w, Y=Y, A=0.0);
    v1, Vadjust, Vnoadjust, polp, pollamb, _, _  = viterFirm(agg, p; maxiter=10000, tol=1e-6, printinfo=false)
    pollamb_dense = makedense(Vadjust, Vnoadjust, p, agg)

    # get joint distribution of prices and shocks
    moments, g0, g = genJointDist(polp, pollamb_dense, p; printinfo=false);


    # get implied aggregate Y
    aggp = 0.0
    for pidx in 1:nquad
        for aidx in 1:nquad
            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]
            density = getDensity(pval, aval, moments, g, g0, p) * wquad[pidx] * wquad[aidx]
            aggp += pval ^ (1.0-p.ϵ) * density
        end
    end
    aggp *= pscale * ascale
    aggp = aggp^(1.0/(1.0-p.ϵ))

    # inteprolate price change function
    itp = Interpolations.interpolate(pollamb, (BSpline(Constant()), BSpline(Constant())))
    eitp = Interpolations.extrapolate(itp, Line())
    pollamb_interp = Interpolations.scale(eitp, pgrid, agrid)
    # get profits to give HH
    # get aggregate fixed cost payments
    # get labour demand
    # integrate
    Ld = 0.0
    F = 0.0
    for pidx = 1:p.nquad
        for aidx = 1:p.nquad
            pval = pgrid_quad[pidx]
            aval = agrid_quad[aidx]
            density = getDensity(pval, aval, moments, g, g0, p) * wquad[pidx] * wquad[aidx]
            Ld += pval^(-p.ϵ) * exp(-aval) * Y * density

            # F is done after realizing shocks
            for epsidx = 1:ngh
                epsval = xgh[epsidx]
                a1val = ρ*aval + σ*epsval
                F += p.κ * pollamb_interp(pval, a1val) * density * wgh[epsidx] # who adjusts in a period
            end

        end
    end
    F *= pscale * ascale
    Ld *= pscale * ascale

    C = Y - F
    w_implied = p.ζ * Ld^(1/p.ν) * C

    # sum of squared errors
    errorY = (aggp - 1.0)^2.0
    errorw = (w - w_implied)^2.0
    error = errorY + errorw
    return error, w, Y, Vadjust, Vnoadjust, polp, pollamb, g, g0, moments, C

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
