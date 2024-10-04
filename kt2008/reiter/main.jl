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

Vadjust1 = repeat(log.(params.kgrid), 1, params.nz)
Vnoadjust1 = zeros(params.nk, params.nz)
V0 = repeat(log.(params.kgrid), 1, params.nz)
kpol = zeros(params.nk, params.nz)

T_adjust!(Vadjust1, kpol, V0, params, agg)
T_noadjust!(Vnoadjust1, V0, params, agg)
V1, xibar_mat = getVout(Vadjust1, Vnoadjust1, params)
error = maximum(abs.(V1 - V0))
V0[:] = V1
@show error
plot(params.kgrid, V1[:, 3])
plot(params.zgrid, V1[5, :])

ki = 1
zi = 1
kval = params.kgrid[ki]
zval = exp(params.zgrid[zi])

function objective(k1)
    out = -agg.P * k1  # only term relevant for k

    Ex = 0.0
    for z1i = 1:params.nz
        Vinterp = interpolate(params.kgrid, V0[:, z1i], BSplineOrder(2))
        Vinterp = extrapolate(Vinterp, Smooth())
        Ex += params.zP[zi, z1i] * Vinterp(k1)
    end

    out += params.beta * 1.0 * Ex

    # maximimze
    out = -1.0 * out

    return out
end
result = optimize(objective, params.kmin, params.kmax, Brent())
@show zval
@show kstar = Optim.minimizer(result)
maxtillnow = -1.0 * Optim.minimum(result)
@show maxtillnow += agg.P*(freduced(exp(params.zgrid[1]), params.kgrid[1], agg, params) + (1.0 - params.delta)*params.kgrid[1])

plot(params.kgrid, objective.(params.kgrid))

plot(params.kgrid, freduced(exp(params.zgrid[1]), params.kgrid, agg, params))


Vadjust, Vnoadjust, Vout, kpol, xibar_mat, error, iter = viter(agg, params)


plot(exp.(params.zgrid), Vadjust[5, :], label="adjust")
plot!(exp.(params.zgrid), Vnoadjust[5, :], label="noadjust")
plot!(exp.(params.zgrid), Vout[5, :], label="overall")