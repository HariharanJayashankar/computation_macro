using MKL
using Plots
using JLD2
using Optim


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
    np_fine = 200, # number of price grids on histogram
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
    phi = 2.2*pflex,
    pgrid = range(plo, phi, length=np),
    pgrid_fine = range(plo, phi, length=np_fine),
    ρ_agg = 0.9
)
p = params()


# equilibrium
# w, Y, Vadjust, Vnoadjust, polp, pollamb, omega, omegahat, C, iter, error = findEquilibrium_ss(p, winit=p.w_flex, Yinit=0.5; deltaY=0.05, deltaw=0.05)
result = optimize(x -> equilibriumResidual(x, p)[1], [1.0, 1.0])
w, Y = Optim.minimizer(result)
error, w, Y, Vadjust, Vnoadjust, polp, pollamb, omega, omegahat, C = equilibriumResidual([w,Y], p)
pdist = sum(omega, dims=2)
heatmap(p.agrid, p.pgrid_fine,omega)
plot(p.pgrid_fine, pdist)
Vout = max.(repeat(Vadjust', p.np, 1), Vnoadjust)

pollamb_mat = zeros(p.np_fine, p.na)
for pidx=1:p.np_fine
    for aidx=1:p.na
        pollamb_mat[pidx, aidx] = pollamb(p.pgrid_fine[pidx], p.agrid[aidx])
    end
end
heatmap(pollamb_mat)

# testing reiter resid
sizedist = p.na * p.np_fine
sizev = p.na * p.np
xss = [
    omega...,
    Vout..., # Vadjust only varies by a
    w,
    p.iss,
    Y,
    C,
    1.0,
    1e-9
]
ηss = zeros(sizev+ 1)
ϵ_ss = zeros(2)

Fout = residequations(xss, xss, ηss, ϵ_ss, p, Y)
@show maximum(abs.(Fout))
findall(abs.(Fout) .≈ maximum(abs.(Fout)))
H1 = FiniteDiff.finite_difference_jacobian(t -> residequations(t, xss,  ηss, ϵ_ss, p, Y), xss)
H2 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, t,  ηss, ϵ_ss, p, Y), xss)
H3 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  t, ϵ_ss, p, Y), ηss)
H4 = FiniteDiff.finite_difference_jacobian(t -> residequations(xss, xss,  xss, t, p, Y), ϵ_ss)


println("Running Gensys")
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = gensysdt(-H2, H1,zeros(size(xss,1)), H4, H3)
@show eu



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

irf_vars = irf[(end-5):end, :]

pw = plot(1:(Tirf), irf_vars[1, :], title="Wage")
pr = plot(1:(Tirf), irf_vars[2, :], title="Interest Rate")
pY = plot(1:(Tirf), irf_vars[3, :], title="Output")
pC = plot(1:(Tirf), irf_vars[4, :], title="Consumption")
pZ = plot(1:(Tirf), irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="Monetary Policy")
pinfl = plot(1:(Tirf), 100.0 * (exp.(irf_vars[6, :]) .- 1.0), title="Inflation (%)")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("reiter2008/tfp_shock.png")

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

irf_vars = irf[(end-5):end, :]
pw = plot(1:Tirf, irf_vars[1, :], title="Wage")
pr = plot(1:Tirf, irf_vars[2, :], title="Interest Rate")
pY = plot(1:Tirf, irf_vars[3, :], title="Output")
pC = plot(1:Tirf, irf_vars[4, :], title="Consumption")
pZ = plot(1:Tirf, irf_vars[5, :], title="TFP")
# pZmon = plot(1:Tirf, irf_vars[6, :], title="TR Shock")
pinfl = plot(1:Tirf, 100.0 * (exp.(irf_vars[6, :]) .- 1.0), title="Inflation (%)")
plot(pw, pr, pY, pC, pZ, pinfl, layout=(2,4), legend=false)
savefig("reiter2008/mon_shock.png")
