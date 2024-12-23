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


param_gen = @with_kw (
    curv = 0.5,
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
    m =  2, # tauchen grid distance
    na = 5, #number of grids in shock
    np = 30, # number of price grids
    npdense = 50, # number of price grids
    γ = 0.05, # learning rte for equilibrium
    # getting shock grid
    shock = tauchen(na, ρ, σ, 0.0, m),
    # shock = rouwenhorst(na, ρ, σ, 0.0),
    aP = shock.p,
    aPstationary = findStationary(aP),
    agrid = shock.state_values,
    flexsoln = flexsoln(ϵ, agrid, aPstationary, ζ, ν),
    pflex = flexsoln[1],
    Pdf_flex = flexsoln[2],
    w_flex = flexsoln[3],
    L_flex = flexsoln[4],
    Y_flex = flexsoln[5],
    plo = 0.5*pflex,
    phi = 1.5*pflex,
    plo_w = 0.01*pflex,
    phi_w =  4.0*pflex,
    pgrid = range(plo^curv, phi^curv, length=np) .^ (1.0/curv),
    pgrid_dense = range(plo^curv, phi^curv, length=npdense) .^ (1.0/curv),
    ρ_agg = 0.9,
    ng = 5,
    # quadrature stuff for winberry
    dampening = 0.1,
    nsimp = 10
)
p = param_gen(ng=1, dampening=1.0)

## does iterateDist work fine?

agg = (w=1.0, Y=1.0, A=0.0);
v1, Vadjust, Vnoadjust, polp, pollamb, iter, err  = viterFirm(agg, p; maxiter=200, tol=1e-6)

pollamb_dense = makedense(Vadjust, Vnoadjust, p, agg)

# get joint distribution of prices and shocks
moments, g0, g, omega1, omega1hat = genJointDist(polp, pollamb_dense, p, printinterval=5);

## ploting policies
plot(p.agrid, polp)

heatmap(pollamb_dense)

## plots of joint distribution
jointdist = zeros(p.npdense, p.na)
for pidx = 1:p.npdense
    for aidx=1:p.na
        pval = p.pgrid_dense[pidx]
        g0_i = g0[aidx]
        g_i = g[:, aidx]
        moments_i = moments[:, aidx]
        
        jointdist[pidx, aidx] = getDensity(pval, moments_i, g_i, g0_i, p)
    end
end

plot(p.pgrid_dense, jointdist[:, 1], label="A lo", title="P dist")
plot!(p.pgrid_dense, jointdist[:, Int(round(end/2))], label="A mid", title="P dist")
plot!(p.pgrid_dense, jointdist[:, end], label="A high", title="P dist")

savefig("ss_dist.pdf")

## Plotting discretized dists
plot(p.pgrid_dense, omega1hat[:, 1], label="A lo", title="P dist")
plot!(p.pgrid_dense, omega1hat[:, 5], label="A mid", title="P dist")
plot!(p.pgrid_dense, omega1hat[:, 10], label="A high", title="P dist")

## Testing parametrization
m0 = rand(p.ng)
m0[1] = 0.12
# m0[2] = 0.0005
gprev = zeros(p.ng)

result = optimize(
    x -> objectiveDensity(x, m0, p),
    (G,x) -> getDensityG!(G, x, m0, p),
    gprev,
    LBFGS(),
    Optim.Options(iterations=1_000_000, show_trace=true)
)
gest = Optim.minimizer(result)
densityOut = Optim.minimum(result)
H = ForwardDiff.hessian(x -> objectiveDensity(x, m0, p), gest)
evals = eigvals(H)
conditionnumber = maximum(evals) / minimum(evals)




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
