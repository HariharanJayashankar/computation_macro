using MKL
using Plots
using JLD2

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
    plo = 0.5*pflex,
    phi = 2.0*pflex,
    pgrid = range(plo, phi, length=np),
    ρ_agg = 0.9,
    ng = 2,
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
    dampening = 0.01
)
p = params()

## Does my gradient accruately give the correct value
m0 = rand(p.ng)
m0 = m0 / 100
m0[1] = log(p.phi) / 4.0

gprev = ones(p.ng) * 0.02
objval = objectiveDensity(gprev, m0, p)
G_fwd = ForwardDiff.gradient(g -> objectiveDensity(g, m0, p), gprev)

G_manual = zeros(p.ng)
getDensityG!(G_manual, gprev, m0, p)
maximum(abs.(G_fwd - G_manual))
# perfect with simpson quad

## Testing dist parametrization


# winberry moments for ng = 2
m0 = rand(p.ng)
m0[1] = log(p.phi) / 4.0
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


## reconstructing moments using the parameters

pmean = simps(x -> log(x) * getDensity(x, m0, gest, 1.0/densityOut, p), p.plo, p.phi, p.nsimp)
# exactly!


