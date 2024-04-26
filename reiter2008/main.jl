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
    iss = 1.0/β - 1.0,
    # otehr parameters (numeric mostly)
    m =  3, # tauchen grid distance
    na = 50, #number of grids in shock
    np = 500, # number of price grids
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
    ρ_agg = 0.9
)
p = params()

vout, vadjust, vnoadjust, polpout, pollambdaout, iter, error = viterFirm((w=p.w_flex, Y=4), p)

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

# equilibrium
w, Y, Vadjust, Vnoadjust, polp, pollamb, omega, omegahat, C, iter, error = findEquilibrium_ss(p, winit=p.w_flex, Yinit=p.Y_flex)


# testing reiter resid
sizedist = p.na * p.np
xss = [
    vec(omegahat)...,
    vec(omega)...,
    vec(Vadjust)...,
    vec(Vnoadjust)...,
    log(w),
    log(p.iss),
    log(Y),
    log(C),
    0.0,
    0.0
]
ηss = zeros(2*sizedist+2)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, p)
H1 = ForwardDiff.jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, p), xss)
# H2 = ForwardDiff.jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, params), xss)
# H3 = ForwardDiff.jacobian(t -> residequations(xss, xss, t,  ϵ_ss, params), η_ss)
# H4 = ForwardDiff.jacobian(t -> residequations(xss, xss, ηss, t,  params), ϵ_ss)
H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, params), xss)
H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, params), xss)
H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, params), ηss)
H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, params), ϵ_ss)
