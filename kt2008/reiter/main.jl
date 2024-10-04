# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")

# testing firm problem
params = paramgen(kmin=0.1, kmax=8.0)
pinit = 2.5
winit = params.phi/pinit
agg = (A=1.0, W=winit, P=pinit)

Vadjust1 = zeros(params.nk, params.nz)
Vnoadjust1 = zeros(params.nk, params.nz)
V0 = zeros(params.nk, params.nz)
V20 = zeros(params.nk, params.nz)
kpol = repeat(params.kgrid, 1, params.nz)

T_adjust_given!(Vadjust1, kpol, V0, V20, params, agg)
T_noadjust!(Vnoadjust1, V0,V20, params, agg)
V0, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)
plot(params.kgrid, V0[:, 3])
plot(params.zgrid, V0[5, :])

ki = 1
zi = 1


Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)


plot(exp.(params.zgrid), Vadjust[5, :], label="adjust")
plot!(exp.(params.zgrid), Vnoadjust[5, :], label="noadjust")
plot!(exp.(params.zgrid), Vout[5, :], label="overall")