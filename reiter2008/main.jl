using  Plots
using JLD2

include("menucost_funcs.jl")
include("gensysdt.jl")

# read jacobians
read_jacob = true

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
Vout = max.(Vadjust, Vnoadjust)


# testing reiter resid
sizedist = p.na * p.np
xss = [
    vec(omega)...,
    vec(Vout)...,
    log(w),
    log(p.iss),
    log(Y),
    log(C),
    1e-9,
    1e-9,
    1e-9
]
# droptol!(xss, 1e-10)
ηss = zeros(sizedist+2)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, p)
if !read_jacob
    H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, p), xss)
    H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, p), xss)
    H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, p), ηss)
    H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, p), ϵ_ss)
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
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = gensysdt(-H1, H2,zeros(size(xss,1)), H3, H4);


println("Making plots...")
# ==  plot IRFS == #
Tirf = 40

# TFP SHOCK
ϵ_tfp_irf = zeros(2, Tirf)
ϵ_tfp_irf[1, 1] = 0.01 # shock value for tfp
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        # everything is deviations from ss 
        irf[:, t] = impact *  ϵ_tfp_irf
    else
        irf[:, t] = G1 * irf[:, t-1]
    end
end

wirf, rirf, Yirf, Cirf, Zirf, Zmonirf  = exp.(irf[(2*sizedist+1):(end-1), :])
infl_irf = irf[end, :]

pw = plot(1:Tirf, wirf, title="Wage")
pr = plot(1:Tirf, rirf, title="Interest Rate")
pY = plot(1:Tirf, Yirf, title="Output")
pC = plot(1:Tirf, Cirf, title="Consumption")
pZ = plot(1:Tirf, Zirf, title="TFP")
pZmon = plot(1:Tirf, Zmonirf, title="Monetary Policy")
pinfl = plot(1:Tirf, inflirf, title="Inflation")
plot(pw, pr, pY, pC, pZ, pZmon, pinfl, layout=(2,4))
savefig("tfp_shock.png")

# monetary policy SHOCK
ϵ_mon_irf = zeros(2, Tirf)
ϵ_mon_irf[2, 1] = 0.01 # shock value for tfp
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        irf[:, t] = impact *  ϵ_mon_irf
    else
        irf[:, t] = G1 * irf[:, t-1]
    end
end

wirf, rirf, Yirf, Cirf, Zirf, Zmonirf  = exp.(irf[(2*sizedist+1):(end-1), :])
infl_irf = irf[end, :]

pw = plot(1:Tirf, wirf, title="Wage")
pr = plot(1:Tirf, rirf, title="Interest Rate")
pY = plot(1:Tirf, Yirf, title="Output")
pC = plot(1:Tirf, Cirf, title="Consumption")
pZ = plot(1:Tirf, Zirf, title="TFP")
pZmon = plot(1:Tirf, Zmonirf, title="Monetary Policy")
pinfl = plot(1:Tirf, inflirf, title="Inflation")
plot(pw, pr, pY, pC, pZ, pZmon, pinfl, layout=(2,4))
savefig("mon_shock.png")
