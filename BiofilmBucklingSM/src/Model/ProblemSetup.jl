function generate_homogeneous_population(T1::TissueMechProperties_t, N, m)
    M = N * m
    k = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    η = ElasticArray{Float64}(zeros(M))

    for i in 1:M
        k[i] = T1.k 
        a[i] = T1.a 
        η[i] = T1.η 
    end

    return BiofilmBucklingSM.HeterogeneousTissueMechProperties_t(k, η, a, T1.restoring_force)
end

function generate_homogeneous_population_LG(T1::TissueMechProperties_LG_t, N, m)
    M = N * m
    k_init = ElasticArray{Float64}(zeros(M))
    k_max = ElasticArray{Float64}(zeros(M))
    λ = ElasticArray{Float64}(zeros(M))
    t_0 = ElasticArray{Float64}(zeros(M))
    a = ElasticArray{Float64}(zeros(M))
    η = ElasticArray{Float64}(zeros(M))

    for i in 1:M
        k_init[i] = T1.k_init 
        k_max[i] = T1.k_max 
        λ[i] = T1.λ 
        t_0[i] = 0.0 
        a[i] = T1.a 
        η[i] = T1.η 
    end

    return BiofilmBucklingSM.HeterogeneousTissueMechProperties_LG_t(k_init, k_max, λ, t_0, η, a, T1.restoring_force)
end


function SetupODEproblem(TissueMech, SimTime, u0, y)
    p = (TissueMech, SimTime, y)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(BiofilmBucklingSM_ODE!,u0,tspan,p), p
end


## to be fixed
function SetupODEproblem_LG(M, TissueMech, SimTime)
    u0 = ElasticArray(collect(LinRange(0,2, M)))
    p = (TissueMech, SimTime)
    tspan = (0.0, SimTime.Tmax)
    return ODEProblem(BiofilmBucklingSM_spring_LG_ODE!,u0,tspan,p), p
end

