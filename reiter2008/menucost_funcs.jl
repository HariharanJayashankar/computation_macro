using QuantEcon, LinearAlgebra, Roots, Parameters

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

    Pdf = (sum(  (exp.(agrid)) .^ (ϵ-1) .* aPstationary))^(1/(ϵ-1)) # price disp
    wout = ( (ϵ-1)/ϵ ) * (1/Pdf)
    pout = sum(((ϵ/(ϵ - 1)) .* wout ./ (exp.(agrid))) .* aPstationary);

    aggPi(L) = L * ( (ϵ/(ϵ-1))^(1-ϵ) - (ϵ/(ϵ-1))^(-ϵ)  ) * wout^(1-ϵ) * Pdf^(ϵ-2)
    consumer_lab(L) = wout / (wout * L - aggPi(L)) - ζ*L^(1/ν)
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
function  viterFirm(agg, params; maxiter=10000, tol=1e-6, printinterval=1000, printinfo=true)


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
        vnoadjust = profit_mat + β * ev0
        vadjust_val = profit_mat .- κ +  β * ev0
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
function  vBackwardFirm(agg, params, v1; maxiter=10000, tol=1e-6, printinterval=1000, printinfo=true)


    @unpack np, na, pgrid, agrid, ϵ, β, aP, κ = params

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


    ev1 = v1 * aP';

    # iterate over choices
    vnoadjust = profit_mat + β * ev1
    vadjust_val = profit_mat .- κ +  β * ev1
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
function Tfunc(omega0hat, polp, pollamb, params)

    omega1 = zeros(params.np, params.na);
    aP = params.aP

    for pidx = 1:params.np
        for aidx = 1:params.na
            
            p1idx = polp[pidx, aidx];

            # non adjusters
            omega1[pidx, aidx] = omega1[pidx, aidx] + (!pollamb[pidx, aidx]) * omega0hat[pidx, aidx];
            
            # adjusters
            omega1[p1idx, aidx] = omega1[p1idx, aidx] + pollamb[pidx, aidx] * omega0hat[pidx, aidx];
        end
    end

    # update shock dist
    omega1hat = zeros(params.np, params.na);

    for pidx = 1:params.np
        for aidx = 1:params.na

            for a0idx = 1:params.na
                omega1hat[pidx, aidx] = omega1hat[pidx, aidx] + omega1[pidx, a0idx] * aP[a0idx, aidx];
            end
        end
    end
    return omega1, omega1hat

end

function genJointDist(polp, pollamb, params; maxiter=1000, tol=1e-6, printinterval=100, printinfo=true)


    omega1hat = ones(params.np, params.na);
    omega1hat = omega1hat ./ (params.np*params.na);
    omega1 = zeros(params.np, params.na);
    error = 10;
    iter = 0;
    
    while (error > tol) && (iter < maxiter)

        omega0hat = omega1hat
        omega1, omega1hat = Tfunc(omega0hat, polp, pollamb, params)
        error = maximum(abs.(omega1hat - omega0hat))
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
function findEquilibrium_ss(p; winit=1, tol=1e-3, max_iter=100, deltaw=0.1,
                            Yinit=1, deltaY=0.1,
                            printinterval=10)

    # gettngg value funcion
    w0 = winit;
    Y0 = Yinit;
    iter = 0;
    error = 10

    # preallocating
    V = zeros(p.np, p.na)
    polp = zeros(p.np, p.na)
    pollamb = zeros(p.np, p.na)
    omega = zeros(p.np, p.na)
    omegahat = zeros(p.np, p.na)

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

    return w0, Y0, V, polp, pollamb,  omega, omegahat, iter, error
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

A lot of it is rewriting the findEquilibrium_ss function
but without iterating till steady state for the distribution
Note r and C are added in compared to steady state since
the euler equation has to hold in equilibrium outside
SS now
Last three are needed to solve todays value function

==#
function residequations(Xl, X, η, ϵ, p)
    
    @unpack np, na = p

    sizedist = np * na

    # unpacking
    omegahat_l = reshape(Xl[1:sizedist], np, na)
    omegahat = reshape(X[1:sizedist], np, na)
    wl, rl, Yl, Cl, Zl, Eyl, Ewl, EMU = Xl[(sizedist+1):end]
    w, r, Y, C, Z, Ey, Ew, Emu = X[(sizedist+1):end]

    # only one shock - tfp shock
    ϵ = ϵ[1]

    # store for reuse
    agrid_orig = p.agrid

    # need to add stochastic discount factors to value function iterations
    # future V
    p.agrid = agrid_orig .* p.ρ_agg .* log(Z)
    V1, Vadjust1, Vnoadjust1, polp1, pollamb1 = viterFirm((w=Ew, Y=Ey), p;
                                        maxiter=10000, tol=1e-6, printinfo=false)
    # current V by backward induction
    p.agrid = agrid_orig .* log(Z)
    V, Vadjust, Vnoadjust, polp, pollamb = vBackwardFirm((w=w, Y=Y), p, V1;
                                        maxiter=10000, tol=1e-6, printinfo=false)
    omega1hat, omega1 = Tfunc(omegahat_l, polp, pollamb, p)
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
            aval = p.agrid[aidx]
            F += p.κ * pollamb[pidx, aidx] .* omegahat[pidx, aidx] # who adjusts in a period
            Ld += pval^(-p.ϵ) * exp(-aval) * Y0 * omega[pidx,aidx]
        end
    end

    cerror = C -  Yimplied - F
    w_implied = p.ζ * Ld^(1/p.ν) * C
    
    euler_error  = (1/C) - (1+r)*p.β*Emu

    # == residual vector == #
    residvector = zeros(np*na + 2, 1)

    # dist error
    omegahat_flat_in = X[1:sizedist]
    omegahat_flat_implied = vec(omega1hat)
    residvector[1:sizedist] = omegahat_flat_in - omegahat_flat_implied

    # other error
    residvector[(sizedist+1):(end-1)] = [Y-Ymplied, w-w_implied]
    residvector[end] = log(Z) - p.ρ_agg * log(Zl) - ϵ

    return residvector
    
end

#==
Linarizeed coefficients
==#
function linearized_coeffs(equations, xss, shocks_sd, fargs)

    shocks_sd = atleast_2d(shocks_sd)
    shocks_ss = zero(shocks_sd)

    H1 = ForwardDiff.jacobian(t -> equations(t, xss, xss, shocks_ss, shocks_sd, fargs...), xss)
    H2 = ForwardDiff.jacobian(t -> equations(xss, t, xss, shocks_ss, shocks_sd, fargs...), xss)
    H3 = ForwardDiff.jacobian(t -> equations(xss, xss, t, shocks_ss, shocks_sd, fargs...), xss)
    H4 = ForwardDiff.jacobian(t -> equations(xss, xss, xss, t, shocks_sd, fargs...), shocks_ss)

    return H1, H2, H3, H4

end