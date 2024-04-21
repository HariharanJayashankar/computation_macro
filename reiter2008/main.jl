using  Plots

include("menucost_funcs.jl")


params = @with_kw (
    β = 0.97^(1/12),
    ζ = 1.0,
    ν = 1.0,
    Π_star = 0.0,
    ϵ = 5.0,
    ρ = 0.8,
    σ = 0.15,
    κ = 0.02,
    iss = 1/β - 1,
    # otehr parameters (numeric mostly)
    m =  3, # tauchen grid distance
    na = 50, #number of grids in shock
    np = 500, # number of price grids
    γ = 0.05, # learning rte for equilibrium
    # getting shock grid
    shock = tauchen(na, ρ, σ, 0, m),
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
    pgrid = range(0.5*pflex, 2*pflex, length=np),
    ρ_agg = 0.9
)
p = params()

@time viterFirm((w=p.w_flex, Y=1), p)
@profview viterFirm((w=p.w_flex, Y=p.Y_flex), p)
vout, vadjust, vnoadjust, polpout, pollambdaout, iter, error = viterFirm((w=p.w_flex, Y=4), p)
vout, vadjust, vnoadjust, polpout, pollambdaout, iter, error = vBackwardFirm((w=p.w_flex, Y=4), p, vout)

heatmap(vout)
heatmap(polpout)
heatmap(pollambdaout)

omegahat, omega = genJointDist(polpout, pollambdaout, p)
heatmap(omegahat)
heatmap(omega)
phatdist = sum(omegahat, dims=2)
plot(p.pgrid, phatdist)
pdist = sum(omega, dims=2)
plot!(p.pgrid, pdist)


w, Y, V, polp, pollamb, omega, omegahat, iter, error = findEquilibrium_ss(p, winit=p.w_flex, Yinit=p.Y_flex)

# 