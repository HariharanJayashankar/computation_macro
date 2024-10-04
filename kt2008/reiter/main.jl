# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")

# testing firm problem
params = paramgen(kmin=0.1, kmax=8.0)
pinit = 2.3
winit = params.phi/pinit
agg = (A=1.0, W=winit, P=pinit)

Vadjust1 = zeros(params.nk, params.nz)
Vnoadjust1 = zeros(params.nk, params.nz)
# V0 = zeros(params.nk, params.nz)
V0 = repeat(log.(params.kgrid), 1, params.nz)
# V20 = zeros(params.nk, params.nz)
V20 = repeat(-1.0./(params.kgrid.^2.0), 1, params.nz)
kpol = repeat(params.kgrid, 1, params.nz)

T_adjust_given!(Vadjust1, kpol, V0, V20, params, agg)
T_noadjust!(Vnoadjust1, V0,V20, params, agg)
V1, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)
V0 = V1
@show Vadjust1
@show Vnoadjust1
# respline up the second derivative
# for zi = 1:params.nz
#     V20[:, zi] = spline(params.kgrid, V0[:, zi], params.nk, 1.0e30, 1.0e30)
# end
plot(params.kgrid, V0[:, 3])
# plot(params.zgrid, V0[5, :])

ki = 1
zi = 1

T_adjust_max!(Vadjust1, kpol, V0, V20, params, agg)
T_noadjust!(Vnoadjust1, V0,V20, params, agg)
V1, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)
V0 = V1
@show Vadjust1
@show Vnoadjust1
# respline up the second derivative
# for zi = 1:params.nz
#     V20[:, zi] = spline(params.kgrid, V0[:, zi], params.nk, 1.0e30, 1.0e30)
# end
plot(params.kgrid, V0[:, 3], label="overall")
plot!(params.kgrid, Vadjust1[:, 3], label="adjust")
plot!(params.kgrid, Vnoadjust1[:, 3], label="noadjust")
# plot(params.zgrid, V0[5, :])


Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)


plot(exp.(params.zgrid), Vadjust[5, :], label="adjust")
plot!(exp.(params.zgrid), Vnoadjust[5, :], label="noadjust")
plot!(exp.(params.zgrid), Vout[5, :], label="overall")