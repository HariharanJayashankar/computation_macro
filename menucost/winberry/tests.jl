using MKL
using Plots
using JLD2

# read jacobians
read_jacob = false

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
    na = 10, #number of grids in shock
    np = 20, # number of price grids
    npdense = 50, # number of price grids
    γ = 0.05, # learning rte for equilibrium
    # getting shock grid
    # shock = tauchen(na, ρ, σ, 0.0, m),
    shock = rouwenhorst(na, ρ, σ, 0.0),
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
    phi = 4.0*pflex,
    pgrid = range(plo^curv, phi^curv, length=np) .^ (1.0/curv),
    pgrid_dense = range(plo^curv, phi^curv, length=npdense) .^ (1.0/curv),
    ρ_agg = 0.9,
    ng = 5,
    # quadrature stuff for winberry
    dampening = 0.1,
    nsimp = 10
)
p = param_gen(ng=2)
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

pmean = simps(x -> x * getDensity(x, m0, gest, 1.0/densityOut, p), p.plo, p.phi, p.nsimp)
@show pmean
@show m0[1]
# exactly!


