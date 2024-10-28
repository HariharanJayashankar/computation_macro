## Set up
using MKL
using Plots
using JLD2
using LinearAlgebra
using ForwardDiff
using Distributions
using Random

include("menucost_funcs.jl")
include("gensysdt.jl")


params = @with_kw (
    β = 0.97^(1/12),
    ζ = 1.0,
    ν = 1.0,
    Π_star = 0.0,
    ϵ = 5.0,
    ρ = 0.8,
    σ = 0.15,
    κ = 0.02,
    iss = 1.0/β - 1.0,
    ϕ_infl = 1.2, # talor rule infl
    ϕ_output = 0.5, # taylor rule output
    # otehr parameters (numeric mostly)
    m =  3, # tauchen grid distance
    na = 10, #number of grids in shock
    np = 20, # number of price grids
    npdense = 50, # number of price grids
    γ = 0.05, # learning rte for equilibrium
    # getting shock grid
    shock = tauchen(na, ρ, σ, 0.0, m),
    aP = shock.p,
    aPstationary = findStationary(aP),
    agrid = shock.state_values,
    flexsoln = flexsoln(ϵ, agrid, aPstationary, ζ, ν),
    pflex = flexsoln[1],
    Pdf_flex = flexsoln[2],
    w_flex = flexsoln[3],
    L_flex = flexsoln[4],
    Y_flex = flexsoln[5],
    plo = 0.1*pflex,
    phi = 2.5*pflex,
    pgrid = range(plo, phi, length=np),
    pgrid_dense = range(plo, phi, length=npdense),
    ρ_agg = 0.9,
    ng = 5,
    nparams = ng > 1 ? 2 + sum([(i+1 for i=2:ng)...]) : 2,
    nquad = 10,
    ngh=3,
    # quadrature stuff for winberry
    quadrature_out = gausslegendre(nquad),
    xquad = quadrature_out[1],
    wquad = quadrature_out[2],
    # for shocks
    gh_quad_out = gausshermite(ngh; normalize=true),
    xgh = gh_quad_out[1],
    wgh = gh_quad_out[2],
    dampening = 0.01,
    nsimp = 10
)
p = params(ng=2)

## Testing
# does simpson2d work well
# https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/math_integration_2D.pdf
# should be 416
f(x,y) = x^2 * y^3
simps2d((x,y) -> f(x,y), 0.0, 2.0, 1.0, 5.0, 6, 6)

# using  gl quadruatrue method
xgrid_quad = scaleUp(p.xquad, 0.0, 2.0)
xscale = (2.0 - 0.0)/2.0
ygrid_quad = scaleUp(p.xquad, 1.0, 5.0)
yscale = (5.0 - 1.0)/2.0
fout = 0.0
for xidx = 1:p.nquad
    for yidx = 1:p.nquad
        xval = xgrid_quad[xidx]
        yval = ygrid_quad[yidx]
        fout += f(xval, yval) * p.wquad[xidx] * p.wquad[yidx]
    end
end

fout *= xscale * yscale
@show fout
# this works!
# simp2d and the gauss legendre quad work.



## Does my gradient accruately give the correct value
m0 = rand(p.nparams)
m0 = m0 / 100
m0[1] = log(p.phi) / 4.0

gprev = ones(p.nparams) * 0.02
objval = objectiveDensity(gprev, m0, p)
G_fwd = ForwardDiff.gradient(g -> objectiveDensity(g, m0, p), gprev)

G_manual = zeros(p.nparams)
getDensityG!(G_manual, gprev, m0, p)
maximum(abs.(G_fwd - G_manual))
# perfect with simpson quad

## Does my integration method integrate over pdfs well?

Random.seed!(123)
d = MvNormal([0.0, 0.0], I)
pdf(d, [0.5, 0.5])
integral = simps2d((x,y) -> pdf(d, [x,y]), -3.0, 3.0, -3.0, 3.0, p.nsimp, p.nsimp)
# not far off from 1

## Testing dist parametrization


# winberry moments for ng = 2
m0 = [
    1.7014
    0.0000
    0.0009
    0.0053
    0.0508]
# m0 = rand(p.nparams)
# m0 = m0 / 100
# m0[1] = log(p.phi) / 4.0
# gprev = ones(p.nparams) * 10.0
# gprev = rand(p.nparams) * -10.0
gprev = zeros(p.nparams)

result = optimize(
    x -> objectiveDensity(x, m0, p),
    (G,x) -> getDensityG!(G, x, m0, p),
    gprev,
    LBFGS(),
    Optim.Options(iterations=1_000_000, show_trace=true)
)
# result = optimize(
#     TwiceDifferentiable(x -> objectiveDensity(x, m0, p), gprev, autodiff=:forward),
#     gprev,
#     Newton(),
#     Optim.Options(iterations=1_000_000, show_trace=true)
# )
gest = Optim.minimizer(result)
densityOut = Optim.minimum(result)
H = ForwardDiff.hessian(x -> objectiveDensity(x, m0, p), gest)
evals = eigvals(H)
conditionnumber = maximum(evals) / minimum(evals)
# problem is that the minimum is exactly at func value 0
# so cant really normalize
# IMplies something is wrong with the method
# dumb fix - putting a high tol on the minimzer stops it from being too bad

## reconstructing moments using the parameters
pmean = 0.0
amean = 0.0
pgrid_quad = scaleUp(p.xquad, p.plo, p.phi)
agrid_quad = scaleUp(p.xquad, minimum(p.agrid), maximum(p.agrid))
pscale = (p.phi - p.plo)/2.0
ascale = (minimum(p.agrid) - maximum(p.agrid))/2.0
for pidx=1:p.nquad
    for aidx=1:p.nquad
        pval = pgrid_quad[pidx]
        aval = agrid_quad[aidx]
        pmean += log(pval) * getDensity(log(pval), aval, m0, gest, 1.0/densityOut, p) * p.wquad[pidx] * p.wquad[aidx]
        amean += aval * getDensity(log(pval), aval, m0, gest, 1.0/densityOut, p) * p.wquad[pidx] * p.wquad[aidx]
    end
end

pmean += pscale * ascale
amean += pscale * ascale
# not even close


## Main run
# equilibrium
result = optimize(
    x -> equilibriumResidual(x, p)[1], 
    [1.0, 1.0],
    Optim.Options(g_tol=1e-4,show_trace=true)
)
w, Y = Optim.minimizer(result)
error, w, Y, Vadjust, Vnoadjust, polp, pollamb, g, g0, moments, C = equilibriumResidual([w,Y], p)
# w0, Y0, Vadjust, Vnoadjust, polp, pollamb,  moments, g, g0, C, iter, error = findEquilibrium_ss(p, winit=p.w_flex, Yinit=p.Y_flex)
Vout = max.(Vadjust, Vnoadjust)

## 
# plotting joint dist with fine grids
pgrid_fine = range(p.plo, p.phi, length=1000)
alo = minimum(p.agrid)
ahi = maximum(p.agrid)
agrid_fine = range(alo, ahi, length=1000)
dist = zeros(1000, 1000)
for i=1:1000
    for j=1:1000
        pval = pgrid_fine[i]
        aval = agrid_fine[j]
        dist[i,j] = getDensity(pval, aval, moments, g, g0, p)
    end
end

heatmap(exp.(agrid_fine), pgrid_fine, dist, xlabel="Shock", ylabel="Price")


## testing reiter resid
sizeval = p.na * p.np
sizedist = p.nparams
xss = [
    g0,
    g...,
    moments...,
    Vadjust...,
    Vnoadjust...,
    w0,
    p.iss,
    Y0,
    C,
    1.0,
    1e-9
]
ηss = zeros(sizeval+1)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, p, Y0)
H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, p, Y0), xss)
H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, p, Y0), xss)
H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, p, Y0), ηss)
H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, p, Y0), ϵ_ss)

println("Running Gensys")
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = gensysdt(-H2, H1,zeros(size(xss,1)), H4, H3)


println("Making plots...")
# ==  plot IRFS == #
Tirf = 20

# TFP SHOCK
ϵ_tfp_irf = zeros(2, Tirf)
ϵ_tfp_irf[1, 1] = p.σ # shock value for tfp
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        # everything is deviations from ss 
        irf[:, t] =  impact * ϵ_tfp_irf[:, 1]
    else
        irf[:, t] = G1 * irf[:, t-1]
    end
end

irf_vars = irf[(2*sizedist+1):end, :]

pw = plot(1:(Tirf), irf_vars[1, :], title="Wage")
pr = plot(1:(Tirf), irf_vars[2, :], title="Interest Rate")
pY = plot(1:(Tirf), irf_vars[3, :], title="Output")
pC = plot(1:(Tirf), irf_vars[4, :], title="Consumption")
pZ = plot(1:(Tirf), irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="Monetary Policy")
pinfl = plot(1:(Tirf), irf_vars[6, :], title="Inflation")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("tfp_shock.png")

# monetary policy SHOCK
ϵ_mon_irf = zeros(2, Tirf)
ϵ_mon_irf[2, 1] = 0.25 # shock value
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        irf[:, t] = impact *  ϵ_mon_irf[:, 1]
    else
        irf[:, t] = G1 * (irf[:, t-1])
    end
end

irf_vars = irf[(2*sizedist+1):end, :]
pw = plot(1:Tirf, irf_vars[1, :], title="Wage")
pr = plot(1:Tirf, irf_vars[2, :], title="Interest Rate")
pY = plot(1:Tirf, irf_vars[3, :], title="Output")
pC = plot(1:Tirf, irf_vars[4, :], title="Consumption")
pZ = plot(1:Tirf, irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="TR Shock")
pinfl = plot(1:Tirf, 100.0 * (exp.(irf_vars[6, :]) .- 1.0), title="Inflation (%)")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("mon_shock.png")
