
using MKL
using Plots
using JLD2
using Optim


include("menucost_funcs.jl")
include("gensysdt.jl")

# read jacobians
read_jacob = false

paramgen = @with_kw (
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
    npdense = 50, # number of price grids on histogram
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
    plo = 0.3*pflex,
    phi = 2.5*pflex,
    pgrid = range(plo, phi, length=np),
    pgrid_dense = range(plo, phi, length=npdense),
    ρ_agg = 0.9
)
params = paramgen()

# manally testing viter

@unpack np, na, pgrid, agrid, ϵ, β, aP, κ, plo, phi = params
agg = (Y=1.0, w=1.0, A=0.0)
Vadjust1 = zeros(na)
Vnoadjust1 = zeros(np, na)
Vadjust0 = zeros(na)
Vnoadjust0 = zeros(np, na)
polp = collect(range(phi, plo, length=na))
pollamb = BitArray(undef, np, na)
v0 = zeros(np, na)

T_adjust_max!(Vadjust1, polp, Vadjust0, Vnoadjust0, params, agg)

v1, Vadjust1, Vnoadjust1, polp, pollamb, iter, error = viterFirm((w=1.0, Y=1.0, A=0.0), params, howarditer=50)

plot(params.pgrid, v1[:, 3])
plot!(params.pgrid, v1[:, 1])
plot!(params.pgrid, v1[:, 5])

plot(params.agrid, Vadjust1)
plot(params.agrid, Vnoadjust1[5, :])
heatmap(params.agrid, params.pgrid, pollamb')
plot(params.agrid, polp)