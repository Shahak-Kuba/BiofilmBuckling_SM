@kwdef struct TissueMechProperties_t
    k::Float64 = 25.0
    η::Float64 = 1.0
    a::Float64 = 10.0
    restoring_force::String = "Hookean"
end

@kwdef struct TissueMechProperties_LG_t
    k_init = 20.0
    k_max = 1.0
    λ::Float64 = 1.0
    η::Float64 = 1.0
    a::Float64 = 10.0
    restoring_force::String = "Hookean"
end

@kwdef mutable struct HeterogeneousTissueMechProperties_t
    k::ElasticArray{Float64}
    η::ElasticArray{Float64}
    a::ElasticArray{Float64}
    restoring_force::String = "Hookean"
end

@kwdef mutable struct HeterogeneousTissueMechProperties_LG_t
    k_init::ElasticArray{Float64}
    k_max::ElasticArray{Float64}
    λ::ElasticArray{Float64}
    t_0::ElasticArray{Float64}
    η::ElasticArray{Float64}
    a::ElasticArray{Float64}
    restoring_force::String = "Hookean"
end

@kwdef struct SimTime_t
    Tmax::Float64 = 20.0 # days
    δt::Float64 = 0.01
    event_δt::Float64 = 0.01
end