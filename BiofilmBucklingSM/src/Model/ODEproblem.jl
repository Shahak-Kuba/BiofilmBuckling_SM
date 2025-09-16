
function BiofilmBucklingSM_ODE!(du, u, p , t)
    TissueMech, SimTime, y = p
    uᵢ₊₁ = circshift(u,(0,-1))
    uᵢ₋₁ = circshift(u,(0,1))

    l,ρ,dv,τ,n = BiofilmBucklingSM.calc_l_ρ_dv_τ_n(uᵢ₊₁, u, uᵢ₋₁, "inward")

    du[:,1] = [0;0]
    du[:,2:end-1] = (1 ./ TissueMech.η[1]) .* ( BiofilmBucklingSM.F(l[2:end-1],TissueMech.k[2:end-1],TissueMech.a[2:end-1],TissueMech.restoring_force).* dv[:,2:end-1]       -      BiofilmBucklingSM.F(circshift(l,(0,1))[2:end-1],TissueMech.k[2:end-1],TissueMech.a[2:end-1],TissueMech.restoring_force).* circshift(dv,(0,1))[:,2:end-1] )
    du[:,end] = [0;0]

    du = BiofilmBucklingSM.check_dr(u, du, y, SimTime.δt)

    # du[:,2:end-1] = (1 ./ TissueMech.η[1]) .* ( BBSM.F(l[2:end-1],TissueMech.k[2:end-1],TissueMech.a[2:end-1],TissueMech.restoring_force).* dv[:,2:end-1]       -      BBSM.F(circshift(l,(0,1))[2:end-1],TissueMech.k[2:end-1],TissueMech.a[2:end-1],TissueMech.restoring_force).* circshift(dv,(0,1))[:,2:end-1] )
end

function BiofilmBucklingSM_spring_LG_ODE!(du, u, p , t)
    TissueMech, SimTime = p
    u_right = circshift(u,-1)
    u_left = circshift(u,1)

    k = BiofilmBucklingSM.k_log_growth.(TissueMech.k_init, TissueMech.k_max, TissueMech.λ, TissueMech.t_0, t)

    du[1] = 0
    du[2:end-1] = (1 ./ TissueMech.η[2:end-1]) .* (Fs(k[2:end-1], TissueMech.a[2:end-1], u[2:end-1], u_right[2:end-1]) .- Fs(k[2:end-1], TissueMech.a[2:end-1], u[2:end-1], u_left[2:end-1]))
    du[end] = -(1 / TissueMech.η[end]) .* (Fs(k[end-1], TissueMech.a[end], u[end], u_left[end]))
end