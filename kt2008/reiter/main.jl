# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")

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
omegass = ssDist(kpoldense, xipoldense, params)

kdist0 = sum(omega0, dims=2)
kdist1 = sum(omega1, dims=2)
kdist2 = sum(omegass, dims=2)
plot(params.kgrid_dense, kdist0, label="T0")
plot!(params.kgrid_dense, kdist1, label="T1")
# updates somewhat reasonably??




