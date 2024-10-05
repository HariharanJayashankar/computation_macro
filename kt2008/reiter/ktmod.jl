using QuantEcon, LinearAlgebra, Roots, Parameters
# using Interpolations
using ForwardDiff
using FiniteDifferences
using FiniteDiff
using SparseArrays, SparsityDetection
using Optim
using BSplineKit


include("utils.jl")
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
    zgrid = shock.state_values,
    nkdense = 50,
    kgrid_dense = exp.(range(log(kmin), log(kmax), nkdense))
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
Output function optimized for labour
==#
function yreduced(z,k,agg,params)

    @unpack alpha, eta = params
    A = agg.A
    W = agg.W

    exponenteta = 1.0 / (1.0 - eta)

    fout = ((eta/W)^(eta*exponenteta)) * (A*z)^(exponenteta)
    fout = fout * (k^(alpha*exponenteta))
    return fout


end

#==
Labour function
==#
function nreduced(z,k,agg,params)

    @unpack alpha, eta = params
    A = agg.A
    W = agg.W

    exponenteta = 1.0 / (1.0 - eta)

    fout = (eta^exponenteta) * (W^(-1.0*exponenteta))
    fout = fout*(z^exponenteta)*(A^exponenteta)*(k^(alpha*exponenteta))
    return fout


end

#==
Bellman adjustment function for a given policy function
==#
function T_adjust_given!(V, kpol, Vprime, V2prime, params, agg)

    @unpack nk, nz, eta, kgrid, zgrid, alpha, zP, beta , delta, kmin, kmax = params

    P = agg.P

    for zi = 1:nz

        zval = exp(zgrid[zi])

        for ki = 1:nk

            kval = kgrid[ki]
            kstar = kpol[ki, zi]

            val = P*(freduced(zval, kval, agg, params) - kstar + (1.0 - delta)*kval)

            Ex = 0.0
            for z1i = 1:nz
                vnext_fun = extrapolate(interpolate(kgrid, Vprime[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(kstar)
            end

            val += beta *  Ex

            V[ki, zi] = val

        end
    end

end

#==
Bellman operator - conditional on capital adjustment

Used as core in value funciton iteration
Can also be used in reiter residual for
backward induction

V is modified - prefereably enter a zero matrix
same shape as Vprime (tomorrows V)
==#
function T_adjust_max!(V, kpol, Vprime, V2prime, params, agg)

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
                vnext_fun = extrapolate(interpolate(kgrid, Vprime[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(k1)
            end

            out += beta *  Ex

            # maximimze
            out = -1.0 * out

            return out
        end

        result = optimize(objective, kmin, kmax, Brent())
        kstar = Optim.minimizer(result)
        kpol[:, zi] .= kstar

        for ki = 1:nk

            kval = kgrid[ki]

            val = P*(freduced(zval, kval, agg, params))
            val = val - P*(kstar - (1.0 - delta) * kval)

            Ex = 0.0
            for z1i = 1:nz
                vnext_fun = extrapolate(interpolate(kgrid, Vprime[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(kstar)
            end

            val += beta * Ex

            V[ki, zi] = val

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
function T_noadjust!(V, Vprime, V2prime, params, agg)

    @unpack nk, nz, eta, kgrid, zgrid, alpha, zP, beta , delta = params

    P = agg.P
    A = agg.A
    W = agg.W


    for ki = 1:nk
        for zi = 1:nz

            kval = kgrid[ki]
            zval = exp(zgrid[zi])
            # k declines next period if we dont adjust
            kval_adj = max((1.0 - delta)*kval, kgrid[1])

            # optimal is a static problem
            val = P * freduced(zval, kval, agg, params)

            Ex = 0.0
            for z1i = 1:nz
                vnext_fun = extrapolate(interpolate(kgrid, Vprime[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(kval_adj)
           end

            val += beta * Ex

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
function viter(agg, params; 
    tol=1e-6, maxiter=1000, howarditer=50,
    printinterval=100, printinfo=true, learningrate=1.0
    )

    @unpack nk, nz, kmin, kgrid = params

    # initialzie value funcitons
    Vadjust1 = zeros(nk, nz)
    Vnoadjust1 = zeros(nk, nz)
    Vold = zeros(nk, nz)
    V2old = zeros(nk, nz)
    xibar_mat = zeros(nk, nz)
    kpol  = repeat(kgrid, 1, nz)
    xipol = 0.5 * ones(nk, nz)
    error = one(tol) + tol
    iter = 0

    while (iter < maxiter) && (error > tol)
        # howard iter
        # iterate with same policy function
        for hiter = 1:howarditer
            T_adjust_given!(Vadjust1, kpol, Vold, V2old, params, agg)
            T_noadjust!(Vnoadjust1, Vold, V2old, params, agg)
            V1, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)
            Vold = learningrate * V1 + (1.0 - learningrate) * Vold

            # respline it up
            for zi = 1:nz
                V2old[:, zi] = spline(kgrid, Vold[:, zi], nk, 1.0e30, 1.0e30)
            end
        end

        # now adjust policy
        T_adjust_max!(Vadjust1, kpol, Vold, V2old, params, agg)
        T_noadjust!(Vnoadjust1, Vold, V2old, params, agg)
        V1, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)

        error = maximum(abs.(V1 - Vold))
        iter += 1
        Vold = V1

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
Convert from coarse functions to fine functions
- Gives us fine policies for kpol and xibar
- Key calculation is in xibar
==#
function makedense(kpol, V, params, agg)

    @unpack nkdense, nz, phi, nk, kgrid, zgrid, delta, zP, beta, kgrid_dense = params
    P = agg.P


    kpoldense = zeros(nkdense, nz)
    xipol_dense = zeros(nkdense, nz)
    
    # placeholder value funcs
    Vadjust = zeros(nk, nz)
    Vnoadjust = zeros(nk, nz)

    for zi = 1:nz
        kstar = kpol[1, zi] # kpol doesnt vary by k state
        kpoldense[:, zi] .= kstar
        zval = exp(zgrid[zi])

        for ki = 1:nkdense

            kval = kgrid_dense[ki]

            # adjusting
            valadjust = P*(freduced(zval, kval, agg, params) - kstar + (1.0 - delta)*kval)

            Ex = 0.0
            for z1i = 1:nz
                vnext_fun = extrapolate(interpolate(kgrid, V[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(kstar)
            end

            valadjust += beta *  Ex


            # not adjusting
            valnoadjust = P*(freduced(zval, kval, agg, params))
            kval_adj = max((1.0-delta)*kval, kgrid[1])

            Ex = 0.0
            for z1i = 1:nz
                vnext_fun = extrapolate(interpolate(kgrid, V[:, z1i], BSplineOrder(4)), Smooth())
                # vnext = splint(kgrid, Vprime[:, z1i], V2prime[:, z1i], nk, kstar)
                Ex += zP[zi, z1i] * vnext_fun(kval_adj)
            end

            valnoadjust += beta *  Ex

            # getting xi
            xival = (valadjust - valnoadjust)/phi
            xipol_dense[ki, zi] = xival


        end

    end

    return kpoldense, xipol_dense


end
#==
Update distribution by one step
given policy function and previous period distribution

Using Young 2010 stochastic simulation update
==#
function updateDist(omega0, kpol_dense, xipol_dense, params)

    @unpack nk, nz, kmin, kgrid, nkdense, 
    kgrid_dense, xibar, kmin, kmax, zP, delta = params

    omega1 = zeros(nkdense, nz)
    for zi = 1:nz

        kstar = kpol_dense[1, zi]

        for ki =1:nkdense
            xival = xipol_dense[ki, zi]
            probadjust = cdfxi(xival, xibar)

            # adjusters weights
            if kstar > kmin && kstar < kmax
                kidxvals = searchsorted(kgrid_dense, kstar)
                kidx_lo = last(kidxvals)
                kidx_hi = kidx_lo + 1
                total_dist = kgrid_dense[kidx_hi] - kgrid_dense[kidx_lo]

                wt_lo = 1.0 - (kstar - kgrid_dense[kidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo
                

            elseif kstar <= kmin
                kidx_lo = 1
                kidx_hi = 2
                wt_lo = 1.0
                wt_hi = 0.0

            elseif kstar >= kmax
                kidx_hi = nkdense
                kidx_lo = kidx_hi - 1
                wt_lo = 0.0
                wt_hi = 1.0
            end

            for z1i = 1:nz
                omega1[kidx_lo, z1i] += probadjust * zP[zi, z1i] * wt_lo * omega0[ki, zi]
                omega1[kidx_hi, z1i] += probadjust * zP[zi, z1i] * wt_hi * omega0[ki, zi]
            end

            # non adjusters
            kval = kgrid_dense[ki]
            kval_adj = max((1.0 - delta)*kval, kmin)
            if kval_adj > kmin && kval_adj < kmax
                kidxvals = searchsorted(kgrid_dense, kval_adj)
                kidx_lo = last(kidxvals)
                kidx_hi = kidx_lo + 1
                total_dist = kgrid_dense[kidx_hi] - kgrid_dense[kidx_lo]

                wt_lo = 1.0 - (kval_adj - kgrid_dense[kidx_lo])/total_dist
                wt_lo = min(1.0, max(0.0, wt_lo))
                wt_hi = 1.0 - wt_lo
                

            elseif kval_adj <= kmin
                kidx_lo = 1
                kidx_hi = 2
                wt_lo = 1.0
                wt_hi = 0.0

            elseif kval_adj >= kmax
                kidx_hi = nkdense
                kidx_lo = kidx_hi - 1
                wt_lo = 0.0
                wt_hi = 1.0
            end

            for z1i = 1:nz
                omega1[kidx_lo, z1i] += (1.0-probadjust) * zP[zi, z1i] * wt_lo * omega0[ki, zi]
                omega1[kidx_hi, z1i] += (1.0-probadjust) * zP[zi, z1i] * wt_hi * omega0[ki, zi]
            end

        end
    end

    omega1 = omega1 ./ sum(omega1)

    return omega1
end

function ssDist(kpoldense, xipoldense, params; 
    tol=1e-6, maxiter=1000, printinfo=false, printinterval=50)

    @unpack nk, nz, kmin, kgrid, nkdense, 
    kgrid_dense, xibar, kmin, kmax, zP, delta = params

    error = one(tol) + tol
    iter = 0
    omega0 = ones(params.nkdense, params.nz)
    omega0 = omega0 ./ sum(omega0)

    while (iter < maxiter) && (error > tol)

        omega1 = updateDist(omega0, kpoldense, xipoldense, params)
        error = maximum(abs.(omega1 - omega0))
        omega0 = omega1
        iter += 1

        if printinfo && (mod(iter, printinterval) == 0)
            println("SS Dist Iterations $iter, error: $error")
        end

    end

    return omega0


end


#==
Residual of steady state
Gives optimal p (which is what we actually solve for)
==#
function ss_equil_resid(p, params; tol=1e-4, maxiter=100, learningrate=0.1)
    

    @unpack nk, nz, kmin, kgrid, nkdense, zgrid,
        kgrid_dense, xibar, kmin, kmax, zP, delta, phi = params

    winit = phi/p
    agg = (A=1.0, W=winit, P=p)

    Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)
    kpoldense, xipoldense = makedense(kpol, Vout, params, agg)
    omegass = ssDist(kpoldense, xipoldense, params, printinfo=true, printinterval=50)

    # get aggregates
    Yagg = 0.0
    Iagg = 0.0
    Kagg = 0.0
    Nagg = 0.0
    for zi = 1:nz
        for ki = 1:nkdense
            zval = exp(zgrid[zi])
            kval = kgrid_dense[ki]

            kstar = kpoldense[ki, zi]
            xistar = xipoldense[ki, zi]
            prob_adjust = cdfxi(xistar, xibar)

            yi = yreduced(zval, kval, agg, params)
            ival = kstar - (1.0 - delta)*kval

            # integrate
            Yagg += omegass[ki, zi]*yi
            Iagg += omegass[ki, zi]*prob_adjust*ival
            Kagg += omegass[ki, zi]*kval
            Nagg += omegass[ki, zi]*(nreduced(zval,kval,agg,params) + expecxi(xistar, xibar))

        end
    end

    Cagg = Yagg - Iagg
    Cerror = 1/p - Cagg

    return Cerror^2, Vadjust, Vnoadjust, Vout, kpol, omegass, Yagg, Iagg, Nagg

end

#==
Residual of equilibrium system
Used for retier
==#
function Fsys(Xl, X, η, ϵ, params)

    @unpack nk, nz, kmin, kgrid, nkdense, zgrid, beta,
        kgrid_dense, xibar, kmin, kmax, zP, delta, phi, rho_A = params

    sizedist = nkdense * nz
    sizev = nk * nz

    # == Unpacking vecotrs == #
    # unpacking "functions" and distr
    Vadj_l = reshape(Xl[1:sizev], nk, nz)
    Vadj = reshape(X[1:sizev], nk, nz)
    Vnoadj_l = reshape(Xl[sizev+1:2*sizev], nk, nz)
    Vnoadj = reshape(X[sizev+1:2*sizev], nk, nz)
    kpol_l = repeat(Xl[(2*sizev+1):(2*sizev + nz)]', nk, 1)
    kpol = repeat(X[(2*sizev+1):(2*sizev + nz)]', nk, 1)
    omega_l = reshape(Xl[(2*sizev + nz+ 1):(2*sizev + nz+ sizedist)], nkdense, nz)
    omega = reshape(Xl[(2*sizev + nz+ 1):(2*sizev + nz+ sizedist)], nkdense, nz)

    # unpacking aggregates
    pl, Yl, Il, Nl, Al = Xl[(2*sizev + nz+ sizedist +1):end]
    p, Y, I, N, A = X[(2*sizev + nz+ sizedist +1):end]

    # == unpack endogenous exepctation errors
    η_vadj = reshape(η[1:sizev], nk, nz)
    η_vnoadj = reshape(η[sizev+1:2*sizev], nk, nz)
    η_kpol =η[2*sizev+1:end]

    # shock
    epaA = ϵ[1]

    # store aggregates needed by the firm
    agg = (P=p, A=Al, W=phi/p)

    # Backwards V iter
    Vadj_l_check = zero(Vadj_l)
    Vnoadj_l_check = zero(Vnoadj_l)
    V1, xibar_mat = getVout(Vadj, Vnoadj, params)
    T_adjust_given!(Vadj_l_check, kpol_l, V1, V1, params, agg)
    T_noadjust!(Vnoadj_l_check, V1, V1, params, agg)
    Vadj_l_check += η_vadj
    Vnoadj_l_check += η_vnoadj

    # policy func error
    kpol_error = zeros(nz)
    for zi = 1:nz
        Ex = 0.0
        kstar = kpol_l[1, zi]
        for z1i = 1:nz
            vprimeinter = Derivative(1) * extrapolate(interpolate(kgrid, V1[:, z1i], BSplineOrder(4)), Smooth())
            Ex +=  zP[zi, z1i] * vprimeinter(kstar)
        end
        kpol_error[zi] = p - beta*Ex + η_kpol[zi]
    end

    # get xipol_l dense by interpolating value functions
    xipol_l_dense = zeros(nkdense, nz)
    kpol_l_dense = zeros(nkdense, nz)
    for ki = 1:nkdense
        for zi = 1:nz
            zval = exp(zgrid[zi])
            kval = kgrid_dense[ki]
            kstar = kpol_l[1, zi]

            # interpolate lagged value functions
            valadjust_interp = extrapolate(interpolate(kgrid, Vadj_l[:, zi], BSplineOrder(4)), Smooth())
            valadjust = valadjust_interp(kval)
            valnoadjust_interp = extrapolate(interpolate(kgrid, Vnoadj_l[:, zi], BSplineOrder(4)), Smooth())
            valnoadjust = valnoadjust_interp(kval)

            xival = (valadjust-valnoadjust)/phi  
            
            xipol_l_dense[ki, zi] = xival
            kpol_l_dense[ki, zi] = kstar

        end
    end

    # update distributin
    omega_check = updateDist(omega_l, kpol_l_dense, xipol_l_dense, params)
    
    # aggregates
    Yagg = 0.0
    Iagg = 0.0
    Kagg = 0.0
    Nagg = 0.0
    for zi = 1:nz
        for ki = 1:nkdense
            zval = exp(zgrid[zi])
            kval = kgrid_dense[ki]

            kstar = kpol_l_dense[ki, zi]
            xistar = xipol_l_dense[ki, zi]
            prob_adjust = cdfxi(xistar, xibar)

            yi = yreduced(zval, kval, agg, params)
            ival = kstar - (1.0 - delta)*kval

            # integrate
            Yagg += omega_l[ki, zi]*yi
            Iagg += omega_l[ki, zi]*prob_adjust*ival
            Kagg += omega_l[ki, zi]*(prob_adjust*kstar + (1.0-prob_adjust)*(1.0-delta)*kval)
            Nagg += omega_l[ki, zi]*(nreduced(zval,kval,agg,params) + expecxi(xistar, xibar))

        end
    end
    Cagg = Yagg - Iagg

    # aggregate A shocks
    shockerror = log(A) - rho_A * log(Al) - epaA

    # == Fill residual == #
    residual = zero(X)

    residual[1:sizev] = vec(Vadj_l) - vec(Vadj_l_check)
    residual[sizev+1:2*sizev] = vec(Vnoadj_l) - vec(Vnoadj_l_check)
    residual[2*sizev+1:2*sizev+nz] = vec(kpol_error)
    residual[2*sizev+nz+1:2*sizev+nz+sizedist] = vec(omega_check - omega)
    residual[2*sizev+nz+sizedist+1:end] = [
        1/p - Cagg,
        Y - Yagg,
        I - Iagg,
        N - Nagg,
        shockerror
    ]


    return residual

    
end