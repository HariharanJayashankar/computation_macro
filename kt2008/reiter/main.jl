# main file for  Terry's version of kt2008

using Pkg
Pkg.activate(".")
using MKL
using Plots
using Optim

include("ktmod.jl")

# testing firm problem
params = paramgen()
pinit = 0.1
winit = params.phi/pinit
agg = (A=1.0, W=winit, P=pinit)

Vadjust, Vnoadjust, Vold, kpol, xibar_mat, error, iter = viter(agg, params)


plot(exp.(params.zgrid), Vadjust[1, :])
plot!(exp.(params.zgrid), Vnoadjust[1, :])

plot(exp.(params.zgrid), Vout[1, :])