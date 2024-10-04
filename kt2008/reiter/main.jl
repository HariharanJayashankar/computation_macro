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

Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)

plot(exp.(params.zgrid), Vadjust[5, :], label="adjust")
plot!(exp.(params.zgrid), Vnoadjust[5, :], label="noadjust")
plot!(exp.(params.zgrid), Vout[5, :], label="overall")

kpoldense, xipoldense = makedense(kpol, Vout, params, agg)


