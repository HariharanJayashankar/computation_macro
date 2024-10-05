# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")
include("gensysdt.jl")

# testing firm problem
params = paramgen(kmin=0.1, kmax=4.0, nk=10, nkdense=50)
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
kdistss = sum(omegass, dims=2)
shockdistss = sum(omegass, dims=1)


plot(params.zgrid, kpol[1, :])

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

pp = plot(1:(Tirf), 100*irf_vars[1, :]./pss, title="Price (%)")
pY = plot(1:(Tirf), 100*irf_vars[2, :]./Yss, title="Output (%)")
pI = plot(1:(Tirf), 100*irf_vars[3, :]./Iss, title="Investment (%)")
pL = plot(1:(Tirf), 100*irf_vars[4, :]./Nss, title="Labour (%)")
pA = plot(1:(Tirf), irf_vars[5, :], title="TFP")
plot(pp, pY, pI, pL, pA, layout=(2,3), legend=false)
savefig("kt2008/reiter/tfp_shock.png")

# check irf of distribution
irf_kdist = zeros(params.nkdense, Tirf)
irf_shockdist = zeros(params.nz, Tirf)
irf_dist = irf[(2*sizev + params.nz+ 1):(2*sizev + params.nz+ sizedist), :]
for t=1:Tirf
    irf_dist_t = reshape(irf_dist[:, t], params.nkdense, params.nz)
    kdist = sum(irf_dist_t, dims=2)
    shockdist = sum(irf_dist_t, dims=1)
    irf_kdist[:, t] = kdistss + kdist
    irf_shockdist[:, t] = shockdistss + shockdist

end

plot(params.kgrid_dense, irf_kdist[:, 1], label = "Kdist = T=1")
plot!(params.kgrid_dense, irf_kdist[:, 5], label = "Kdist = T=5")
plot!(params.kgrid_dense, irf_kdist[:, 10], label = "Kdist = T=10")
plot!(params.kgrid_dense, irf_kdist[:, Tirf], label = "Kdist = T=final")
plot!(params.kgrid_dense, kdistss, label = "Kdist ss")
# kdist moves DOWN on impact. Not great

plot(params.zgrid, irf_shockdist[:, 1], label = "Shockdist = T=1")
plot!(params.zgrid, irf_shockdist[:, 5], label = "Shockdist = T=5")
plot!(params.zgrid, irf_shockdist[:, 10], label = "Shockdist = T=10")
plot!(params.zgrid, irf_shockdist[:, Tirf], label = "Shockdist = T=final")
plot!(params.zgrid, shockdistss', label = "Shockdist ss")

# lets see how the policy function evolves
irf_polk = repeat(kpol[1, :], 1, 20) + irf_dist[(2*sizev+1):(2*sizev + params.nz), :]
plot(params.zgrid, irf_polk[:, 1], label = "Policy t= 1")
plot!(params.zgrid, irf_polk[:, 5], label = "Policy t= 5")
plot!(params.zgrid, irf_polk[:, 10], label = "Policy t= 10")
plot!(params.zgrid, irf_polk[:, Tirf], label = "Policy t= Final")
plot!(params.zgrid, kpol[1, :], label = "Policy ss")
# policy function doesnt really move over time

# lets check the value functions and adjustment probabilities
# value funcs across z for the middle z value
irf_vaval = zeros(params.nk, Tirf)
irf_vnaval = zeros(params.nk, Tirf)
irf_xi = zeros(params.nk, Tirf)
irf_vaval_sub = irf[1:sizev, :]
irf_vnaval_sub = irf[sizev+1:2*sizev, :]
for t = 1:Tirf
    vaval_t = reshape(irf_vaval_sub[:, t], params.nk, params.nz)
    vaval_t_z = vaval_t[:, 3]
    vnaval_t = reshape(irf_vnaval_sub[:, t], params.nk, params.nz)
    vnaval_t_z = vnaval_t[:, 3]
    xi = (vaval_t_z - vnaval_t_z) ./ params.phi

    irf_vaval[:, t] = vaval_t_z
    irf_vnaval[:, t] = vnaval_t_z
    irf_xi[:, t] = xi
end

plot(params.kgrid, Vadjust[:, 3] + irf_vaval[:, 1], label = "Va t=1")
plot!(params.kgrid,Vadjust[:, 3] +  irf_vaval[:, 5], label = "Va t=5")
plot!(params.kgrid, Vadjust[:, 3] + irf_vaval[:, 10], label = "Va t=10")
plot!(params.kgrid, Vadjust[:, 3], label = "Va ss")

plot(params.kgrid, irf_xi[:, 1], label = "Va t=1")
plot!(params.kgrid, irf_xi[:, 5], label = "Va t=5")
plot!(params.kgrid, irf_xi[:, 10], label = "Va t=10")