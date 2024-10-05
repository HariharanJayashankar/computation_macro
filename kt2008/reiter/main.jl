# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")
include("gensysdt.jl")

# testing firm problem
params = paramgen(kmin=0.1, kmax=8.0)
pinit = 2.3 # a narrrow band of prices only works fori nterior solutions
winit = params.phi/pinit
agg = (A=1.0, W=winit, P=pinit)

Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)

plot(exp.(params.zgrid), Vadjust[5, :], label="adjust")
plot!(exp.(params.zgrid), Vnoadjust[5, :], label="noadjust")
plot!(exp.(params.zgrid), Vout[5, :], label="overall")

kpoldense, xipoldense = makedense(kpol, Vout, params, agg)

# testing distribtuino update
omega0 = ones(params.nkdense, params.nz)
omega0 = omega0 ./ sum(omega0)
omega1 = updateDist(omega0, kpoldense, xipoldense, params)
omegass = ssDist(kpoldense, xipoldense, params, printinfo=true, printinterval=1)

kdist0 = sum(omega0, dims=2)
kdist1 = sum(omega1, dims=2)
kdist2 = sum(omegass, dims=2)
plot(params.kgrid_dense, kdist0, label="T0")
plot!(params.kgrid_dense, kdist1, label="T1")
plot!(params.kgrid_dense, kdist2, label="ss")
# updates somewhat reasonably??

result = optimize(x -> ss_equil_resid(x, params)[1], 2.15, 2.3, GoldenSection())
pss = Optim.minimizer(result)

error, Vadjust, Vnoadjust, Vout, kpol, omegass, Yss, Iss, Nss = ss_equil_resid(pss, params)

sizev = params.nk * params.nz
sizedist = params.nkdense * params.nz
xss = [
    Vadjust...,
    Vnoadjust...,
    kpol[1, :]...,
    omegass...,
    pss,
    Yss,
    Iss,
    Nss,
    1.0
]

ηss = zeros(2*sizev + params.nz)
Fout = Fsys(xss, xss, ηss, [0.0], params)
maximum(abs.(Fout))
findall(abs.(Fout) .≈ maximum(abs.(Fout)))

H1 = FiniteDiff.finite_difference_jacobian(t -> Fsys(t, xss,  ηss, [0.0], params), xss)
H2 = FiniteDiff.finite_difference_jacobian(t -> Fsys(xss, t,  ηss, [0.0], params), xss)
H3 = FiniteDiff.finite_difference_jacobian(t -> Fsys(xss, xss,  t, [0.0], params), ηss)
H4 = FiniteDiff.finite_difference_jacobian(t -> Fsys(xss, xss,  xss, t, params), [0.0])


println("Running Gensys")
G1, Const, impact, fmat, fwt, ywt, gev, eu, loose = gensysdt(-H2, H1,zeros(size(xss,1)), H4, H3)
@show eu

println("Making plots...")
# ==  plot IRFS == #
Tirf = 20

# TFP SHOCK
ϵ_tfp_irf = zeros(1, Tirf)
ϵ_tfp_irf[1, 1] = params.sigma_A # shock value for tfp
irf = zeros(size(xss, 1), Tirf)
for t=1:Tirf
    if t == 1
        # everything is deviations from ss 
        irf[:, t] =  impact * ϵ_tfp_irf[:, 1]
    else
        irf[:, t] = G1 * irf[:, t-1]
    end
end

irf_vars = irf[(end-4):end, :]

pp = plot(1:(Tirf), irf_vars[1, :], title="Price")
pY = plot(1:(Tirf), irf_vars[2, :], title="Output")
pI = plot(1:(Tirf), irf_vars[3, :], title="Investment")
pL = plot(1:(Tirf), irf_vars[4, :], title="Labour")
pA = plot(1:(Tirf), irf_vars[5, :], title="TFP")
plot(pp, pY, pI, pL, pA, layout=(2,3), legend=false
savefig("kt2008/reiter/tfp_shock.png")