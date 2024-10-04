using QuantEcon, LinearAlgebra, Roots, Parameters
# using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using SparseArrays, SparsityDetection
using Optim
using BSplineKit

#==
Get steady state level of capital from neoclassical model
used to calculate grid space
==#
function kss(β, α, δ)
    k = (1/(β*α) + (δ-1)/(α))^(1/(α-1))
    return k
end

# define parameters
# default model params taken from Terry 2017
paramgen = @with_kw (
    alpha = 0.256,
    beta = 0.977,
    delta = 0.069,
    sigma_A = 0.014,
    sigma_z = 0.022,
    eta = 0.640,
    phi = 2.4,
    rho_A = 0.859,
    rho_z = 0.859,
    xibar = 0.0083,
    kmin = 0.2 * kss(beta, alpha, delta),
    kmax = 3.0 * kss(beta, alpha, delta),
    nk = 20,
    kgrid = exp.(range(log(kmin), log(kmax), nk)),
    m =  3, # tauchen grid distance
    nz = 10, #number of grids in shock,
    shock = tauchen(nz, rho_z, sigma_z, 0.0, m),
    zP = shock.p,
    zgrid = shock.state_values
)

# distribution of xi (menu cost)
function cdfxi(xi, xibar)
    if (xi<=0.0)
        out = 0.0
    elseif (xi > 0.0 && xi<=xibar) 
        out = xi/xibar
    elseif (xi>xibar)
        out = 1.0
    end
    return out
end

function expecxi(xi, xibar)
    if xi<=0.0
        out = 0.0
    elseif (xi<=xibar && xi>=0.0)
        out = ( xi ^ 2.0 ) / (2.0 * xibar)
    elseif (xi>xibar)
        out = ( xibar ^ 2.0 ) / (2.0 * xibar)
    end
    return out
end

# production function
function fprod(z,k,n,agg,params)

    @unpack alpha, eta = params

    A = agg.A
    fout = A*z*(k^alpha)*(n^(eta))

    return fout

end

#==
Production function optimized for labour (and its cost)
==#
function freduced(z,k,agg,params)

    @unpack alpha, eta = params
    A = agg.A
    W = agg.W

    exponenteta = 1.0 / (1.0 - eta)

    fout = ( eta ^ (eta * exponenteta) ) * (1.0 - eta) * ( W ^ ( -1.0 * eta * exponenteta ) )
    fout = fout * (z ^ exponenteta) * (A ^ exponenteta) * ( k .^ (alpha * exponenteta) )
    return fout


end

#==
Bellman operator - conditional on capital adjustment

Used as core in value funciton iteration
Can also be used in reiter residual for
backward induction

V is modified - prefereably enter a zero matrix
same shape as Vprime (tomorrows V)
==#
function T_adjust!(V, kpol, Vprime, params, agg; sdf=1.0)

    @unpack nk, nz, eta, kgrid, zgrid, alpha, zP, beta , delta, kmin, kmax = params

    P = agg.P
    A = agg.A
    W = agg.W

    # ev1 = Vprime * zP'

    for zi = 1:nz

        zval = exp(zgrid[zi])
        # using brent to optimize
        # needs Vprime to be an interpolated value function
        # independent of zval

        function objective(k1)
            out = -P * k1  # only term relevant for k

            Ex = 0.0
            for z1i = 1:nz
                Vinterp = interpolate(kgrid, Vprime[:, z1i], BSplineOrder(2))
                Vinterp = extrapolate(Vinterp, Smooth())
                Ex += zP[zi, z1i] * Vinterp(k1)
            end

            out += beta * sdf * Ex

            # maximimze
            out = -1.0 * out

            return out
        end

        result = optimize(objective, kmin, kmax, Brent())
        @show zval
        @show kstar = Optim.minimizer(result)
        maxtillnow = -1.0 * Optim.minimum(result)

        for ki = 1:nk

            kval = kgrid[ki]

            # maxtillnow = -maxtillnow # revert this
            @show maxtillnow += P*(freduced(zval, kval, agg, params) + (1.0 - delta)*kval)
            
            # # max over k
            # maxtillnow = -100000.0
            # kstar = 0
            # for k1i = 1:nk
            #     k1val = kgrid[k1i]

            #     val = P * (freduced(zval, kval, agg, params) - k1val + (1.0 - delta)*kval)
            #     val += beta * sdf * ev1[k1i, zi]

            #     if val > maxtillnow
            #         kstar = k1i
            #         maxtillnow = val
            #     end
            # end

            V[ki, zi] = maxtillnow
            kpol[ki, zi] = kstar

        end
    end

end
#==
Bellman operator - no adjustment

Used as core in value funciton iteration
Can also be used in reiter residual for
backward induction

V is modified - prefereably enter a zero matrix
same shape as Vprime (tomorrows V)
==#
function T_noadjust!(V, Vprime, params, agg; sdf=1.0)

    @unpack nk, nz, eta, kgrid, zgrid, alpha, zP, beta , delta = params

    P = agg.P
    A = agg.A
    W = agg.W

    # k declines next period if we dont adjust
    v1_adj = zeros(eltype(Vprime), nk, nz)
    for zi = 1:nz
        Vinterp = interpolate(kgrid, Vprime[:, zi], BSplineOrder(2))
        Vinterp = extrapolate(Vinterp, Smooth())
        for ki = 1:nk
            kval = kgrid[ki]
            zval = exp(zgrid[zi])
            kval_adj = max((1.0 - delta)*kval, kgrid[1])

            v1_adj[ki, zi] = Vinterp(kval_adj)
        end
    end

    ev1 = v1_adj * zP'

    for ki = 1:nk
        for zi = 1:nz

            kval = kgrid[ki]
            zval = exp(zgrid[zi])

            # optimal is a static problem

            val = P * freduced(zval, kval, agg, params)
            val += beta * sdf * ev1[ki, zi]

            V[ki, zi] = val
        end
    end

end

#==
Get overall value funciton given ones for adjusting and not adjusting
==#
function getVout(Vadjust, Vnoadjust, params)


    @unpack nk, nz, eta, kgrid, zgrid, alpha, zP, beta , delta, phi, xibar = params

    xibar_mat = zeros(nk, nz)
    Vout = zeros(nk, nz)

    for ki = 1:nk
        for zi = 1:nz
            xibartmp = (Vadjust[ki, zi] - Vnoadjust[ki, zi])/phi
            prob_adjust = cdfxi(xibartmp, xibar)

            Vout[ki, zi] = -phi * expecxi(xibartmp, xibar) + 
                prob_adjust * Vadjust[ki, zi] +
                (1.0 - prob_adjust) * Vnoadjust[ki, zi]
            xibar_mat[ki, zi] = xibartmp
        end
    end


    return Vout, xibar_mat

end


#== 
Value function iteration=
==#
function viter(agg, params; tol=1e-6, maxiter=1000, printinterval=100, printinfo=true)

    @unpack nk, nz = params

    # initialzie value funcitons
    Vadjust1 = zeros(nk, nz)
    Vnoadjust1 = zeros(nk, nz)
    Vold = zeros(nk, nz)
    xibar_mat = zeros(nk, nz)
    kpol  = zeros(nk, nz)
    error = one(tol) + tol
    iter = 0

    while (iter < maxiter) && (error > tol)

        Vadjust0 = deepcopy(Vadjust1)
        Vnoadjust0 = deepcopy(Vnoadjust1)
        Vold, xibar_mat = getVout(Vadjust0, Vnoadjust0, params)

        T_adjust!(Vadjust1, kpol, Vold, params, agg)
        T_noadjust!(Vnoadjust1, Vold, params, agg)

        erroradjust = maximum(abs.(Vadjust1 - Vadjust0))
        errornoadjust = maximum(abs.(Vnoadjust1 - Vnoadjust0))
        error = max.(erroradjust, errornoadjust)
        iter += 1

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

    return Vadjust1, Vnoadjust1, Vold, kpol, xibar_mat, error, iter

end

#==
Update distribution by one step
given policy function and previous period distribution
==#
function updateDist(omega0, polk, pollamb)
    return 0
end


#==
Residual of steady state
Gives optimal p (which is what we actually solve for)
==#
function ss_equil_resid(params; tol=1e-4, maxiter=100, learningrate=0.1)


    return 0

end