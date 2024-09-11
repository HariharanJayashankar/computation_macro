using MKL
using Plots
using JLD2
using LinearAlgebra

include("menucost_funcs.jl")
include("gensysdt.jl")

# read jacobians
read_jacob = false

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
    np = 200, # number of price grids
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
    pss = log.(pflex),
    plo = 0.1*pflex,
    phi = 2.0*pflex,
    pgrid = range(plo, phi, length=np),
    ρ_agg = 0.9,
    ng = 5,
    nparams = 2 + sum([(i+1 for i=2:ng)...]),
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
    nsimp=20, # number of points for simpson quadrature
    dampening=0.1 # learning rate
)
p = params()

# does simpson2d work well
# https://d-arora.github.io/Doing-Physics-With-Matlab/mpDocs/math_integration_2D.pdf
# should be 416
f(x,y) = x^2 * y^3
simps2d(f, 0.0, 2.0, 1.0, 5.0, 6, 6)


# manually testing winbrry
# m0 = rand(p.nparams)
# taken from equilibrium with histogram
m0 = [0.9060857580754306, -9.540979117872439e-17, 0.058261844178690333, -0.06243177533634714, 0.06911418673639061, 0.009348699626453561, -0.006537415278756326, 0.0035201241401555147, 4.515394998958746e-17, 0.011542038186751393, -0.011316905912208635, 0.011587607563491472, -0.01238504694776647, 0.013798719971406092, 0.00500138025710313, -0.004023892704450641, 0.003127755735850773, -0.0022234913216252665, 0.001217471208940386, 1.563162908210859e-17]
gprev = rand(p.nparams)
gprev = zeros(p.nparams)
objectiveDensity(gprev, m0, p)
getDensity(1.0, 1.0, m0, gprev, 1.0, p)

# is my gradient correct?
using ForwardDiff
fwd = ForwardDiff.gradient(x -> objectiveDensity(x, m0, p), gprev)
G = zero(m0)
getDensityG!(G, gprev, m0, p)

maximum(abs.(G - fwd))
# seems right


gprev =  ones(p.nparams)
gprev = zeros(p.nparams)
# result = optimize(
#     x -> objectiveDensity(x, m0, p),
#     gprev,
#     Optim.Options(1e-3, g_tol=1e-3, iterations=100_000)
#  )
# result = optimize(
#     x -> objectiveDensity(x, m0, p),
#     gprev,
#     LBFGS(),
#     Optim.Options(x_tol=1e-6, f_tol=1e-6, g_tol=1e-6, iterations=100_000, show_trace=true)
# )
# # Optim.minimizer(result)
result = optimize(
    x -> objectiveDensity(x, m0, p),
    (G,x) -> getDensityG!(G, x, m0, p),
    gprev,
    LBFGS(),
    Optim.Options(x_tol = 1e-3, f_tol = 1e-3, g_tol = 1e-3, iterations=100_000, show_trace=true)
)
@show Optim.converged(result)
gest = Optim.minimizer(result)
# objectiveDensity(gest, m0, p)


# identifiability
# condition number (absolutely) should'nt be bigger than 1
H = ForwardDiff.hessian(x -> objectiveDensity(x, m0, p), gest)
evals = eigvals(H)
conditionnumber = maximum(evals) / minimum(evals)
# problem seems very badly defined!

densityOut = Optim.minimum(result)
g0 = 1.0 / densityOut

# see if this density gives back similar moments
alo = minimum(p.agrid)
ahi = maximum(p.agrid)
meanp = simps2d((price,a) -> price * getDensity(price, a, m0, gest, g0, p), p.plo, p.phi, alo, ahi, p.nsimp, p.nsimp)
meana = simps2d((price,a) -> a * getDensity(price, a, m0, gest, g0, p), p.plo, p.phi, alo, ahi, p.nsimp, p.nsimp)
@show abs(meanp - m0[1])
@show abs(meana - m0[2])


# equilibrium
w0, Y0, Vadjust, Vnoadjust, polp, pollamb,  moments, gmat, g0, C, iter, error = findEquilibrium_ss(p, winit=p.w_flex, Yinit=p.Y_flex)
Vout = max.(Vadjust, Vnoadjust)

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
        dist[i,j] = getDensity(pval, aval, moments, gmat, g0, params)
    end
end

heatmap(dist)


# testing reiter resid
sizeval = p.na * p.np
sizedist = p.ng * (p.ng + 1)
xss = [
    g0,
    gmat...,
    moments...,
    Vout...,
    w,
    p.iss,
    Y,
    C,
    1.0,
    1e-9
]
ηss = zeros(1*sizedist+1)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, p, Y)
if !read_jacob
    H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, p, Y), xss)
    H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, p, Y), xss)
    H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, p, Y), ηss)
    H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, p, Y), ϵ_ss)
    jldsave("solnmatrices.jld2"; H1, H2, H3, H4)
else
    println("Reading Jacobians...")
    H1, H2, H3, H4 = load("solnmatrices.jld2", H1, H2, H3, H4)
    H1 = collect(H1)
    H2 = collect(H2)
    H3 = collect(H3)
    H4 = collect(H4)
end

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
