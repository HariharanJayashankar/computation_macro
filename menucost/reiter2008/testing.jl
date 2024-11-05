## Setup
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
    pss = pflex,
    plo = 0.1*pflex,
    phi = 5.0*pflex,
    pgrid = exp.(range(log(plo), log(phi), length=np)),
    pgrid_dense = range(plo, phi, length=npdense),
    ρ_agg = 0.9
)
params = paramgen()

## Testing
# manally testing viter
@time viterFirm((w=1.0, Y=1.0, A=0.0), params, howarditer=50)
v1, Vadjust1, Vnoadjust1, polp, pollamb, iter, error = viterFirm((w=1.0, Y=1.0, A=0.0), params, howarditer=50)

plot(params.pgrid, v1[:, 3])
plot(params.pgrid, v1[:, 1])
plot(params.pgrid, v1[:, 5])

plot(params.agrid, Vadjust1)
plot(params.agrid, Vnoadjust1[5, :])
heatmap(pollamb')
plot(params.agrid, polp)

pollamb_dense = makedense(Vadjust1, Vnoadjust1, params, (w=1.0, Y=1.0, A=0.0))
omegass, omegahatss = genJointDist(polp, pollamb_dense, params)


plot(params.pgrid_dense, omegahatss[:, 1], label="A lo", title="P dist")
plot!(params.pgrid_dense, omegahatss[:, 3], label="A mid", title="P dist")
plot!(params.pgrid_dense, omegahatss[:, 5], label="A high", title="P dist")

savefig("ss_dist.pdf")
## SS Equilibrium

ssresult = optimize(
    x -> equilibriumResidual(x, params)[1],
    [1.0, 1.0],
    Optim.Options(
        show_trace=true,
        g_tol = 1e-4
    )
)
w, Y = Optim.minimizer(ssresult)
error, w, Y, Vadjust, Vnoadjust, polp, pollamb, omega, omegahat, C = equilibriumResidual(
    [w, Y], params
)

## Linearizing and reiter
sizedist = params.npdense * params.na
sizev = params.np * params.na
xss = [
    omega...,
    Vadjust...,
    Vnoadjust...,
    polp...,
    w,
    params.iss,
    Y,
    C,
    0.0,
    0.0
]


ηss = zeros(2*params.na+sizev+1)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, params, Y)
@show maximum(abs.(Fout))
findall(abs.(Fout) .≈ maximum(abs.(Fout)))
findall(abs.(Fout) .> 1e-2)
H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, params, Y), xss)
H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, params, Y), xss)
H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, params, Y), ηss)
H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, params, Y), ϵ_ss)


println("Running Gensys")
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = gensysdt(-H2, H1,zeros(size(xss,1)), H4, H3)
@show eu

println("Making plots...")
# ==  plot IRFS == #
Tirf = 20

# TFP SHOCK
ϵ_tfp_irf = zeros(2, Tirf)
ϵ_tfp_irf[1, 1] = params.σ # shock value for tfp
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        # everything is deviations from ss 
        irf[:, t] =  impact * ϵ_tfp_irf[:, 1]
    else
        irf[:, t] = G1 * irf[:, t-1]
    end
end

irf_vars = irf[(end-5):end, :]

pw = plot(1:(Tirf), 100*irf_vars[1, :], title="Wage")
pr = plot(1:(Tirf), 100*irf_vars[2, :], title="Interest Rate")
pY = plot(1:(Tirf), 100*irf_vars[3, :], title="Output")
pC = plot(1:(Tirf), 100*irf_vars[4, :], title="Consumption")
pZ = plot(1:(Tirf), irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="Monetary Policy")
pinfl = plot(1:(Tirf), 100.0 * irf_vars[6, :], title="Inflation")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("tfp_shock.png")

# monetary policy SHOCK
ϵ_mon_irf = zeros(2, Tirf)
ϵ_mon_irf[2, 1] = 0.25/100.0 # shock value
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        irf[:, t] = impact *  ϵ_mon_irf[:, 1]
    else
        irf[:, t] = G1 * (irf[:, t-1])
    end
end

irf_vars = irf[(end-5):end, :]
pw = plot(1:Tirf, 100*irf_vars[1, :], title="Wage")
pr = plot(1:Tirf, 100*irf_vars[2, :], title="Interest Rate")
pY = plot(1:Tirf, 100*irf_vars[3, :], title="Output")
pC = plot(1:Tirf, 100*irf_vars[4, :], title="Consumption")
pZ = plot(1:Tirf, irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="TR Shock")
pinfl = plot(1:Tirf, 100.0 * irf_vars[6, :], title="Inflation (%)")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("mon_shock.png")


