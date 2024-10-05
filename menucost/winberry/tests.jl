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


# manually testing winbrry
# m0 = ones(p.nparams)
# gprev = -0.02 * ones(p.nparams)
# gprev = zeros(2 + (p.ng-1)*(p.ng + 1))
m0 = rand(p.nparams)
gprev = -1e-2 * ones(p.nparams)
objectiveDensity(gprev, m0, p)

# is my gradient correct?
using ForwardDiff
fwd = ForwardDiff.gradient(x -> objectiveDensity(x, m0, p), gprev)
G = zero(m0)
getDensityG!(G, gprev, m0, p)

maximum(abs.(G - fwd))
# seems right

# does my integration look righT?
#test for int_-2^2 x^2 dx = 16/3
xlo = -2.0
xhi = 2.0
xgrid_scaled = scaleUp(p.xquad, xlo, xhi)
scaling_factor = (xhi - xlo)/2.0
f(x) = x.^2
(p.wquad ⋅ f.(xgrid_scaled)) * scaling_factor
16/3
# my understandding seems right!



result = optimize(
    x -> objectiveDensity(x, m0, p),
    gprev,
    Optim.Options(x_tol=1e-8, f_tol=1e-8, g_tol=1e-8, iterations=100_000)
 )
result = optimize(
    x -> objectiveDensity(x, m0, p),
    gprev,
    LBFGS(),
    Optim.Options(x_tol=1e-6, f_tol=1e-6, g_tol=1e-6, iterations=100_000, show_trace=true)
)
Optim.minimizer(result)
result = optimize(
    x -> objectiveDensity(x, m0, p),
    (G,x) -> getDensityG!(G, x, m0, p),
    gprev,
    BFGS(),
    Optim.Options(x_tol=1e-6, f_tol=1e-6, g_tol=1e-6, iterations=100_000, show_trace=true)
)
@show Optim.converged(result)
gest = Optim.minimizer(result)
densityOut = Optim.minimum(result)
# g0 = 1.0 / densityOut
g0 = 1.0

# see if this density gives back similar moments
alo = minimum(p.agrid)
ahi = maximum(p.agrid)
pgrid_quad = scaleUp(p.xquad, p.plo, p.phi)
agrid_quad = scaleUp(p.xquad, alo, ahi)
meanp = 0.0
meana = 0.0
for pidx = 1:p.nquad
    for aidx = 1:p.nquad
        pval = pgrid_quad[pidx]
        aval = agrid_quad[aidx]

        density = getDensity(pval, aval, m0, gest, g0, p) * p.wquad[pidx] * p.wquad[aidx]
        meanp += density * pval
        meana += density * aval
    end
end
@show abs(meanp - m0[1])
@show abs(meana - m0[2])
