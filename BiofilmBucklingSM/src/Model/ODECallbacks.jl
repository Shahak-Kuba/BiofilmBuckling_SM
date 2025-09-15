using ElasticArrays

# Modifier functions for insert, push and deleteat for ElasticArrays and ElasticVectors

function Base.insert!(A::ElasticArray, col, val)
    data = A.data
    insert!(data, col, val)
    return A
end

function Base.insert!(A::ElasticVector, col, val)
    data = A.data
    insert!(data, col, val)
    return A
end


function Base.push!(A::ElasticVector, val::String)
    data = A.data
    push!(data, val)
    return A
end


function Base.deleteat!(A::ElasticVector, col)
    data = A.data
    deleteat!(data,col)
    return A
end

function tissue_dep_affect!(integrator)
    u = integrator.u
    TissueMech, SimTime = integrator.p
    #curr_t = integrator.t
    #println("inserting a new point at t = ", integrator.t)
    idx = size(u, 1)
    new_node_ratio = 1.0
    insert!(u, idx, u[end-1] + new_node_ratio *((u[end] - u[end-1])))
    insert!(TissueMech.k, idx, TissueMech.k[idx])
    insert!(TissueMech.η, idx, TissueMech.η[idx])
    insert!(TissueMech.a, idx, TissueMech.a[idx])
    resize!(integrator, size(integrator.u, 1))
    nothing
end

function tissue_dep_LG_affect!(integrator)
    u = integrator.u
    TissueMech, SimTime = integrator.p
    #curr_t = integrator.t
    #println("inserting a new point at t = ", integrator.t)
    idx = size(u, 1)
    new_node_ratio = 0.95
    insert!(u, idx, u[end-1] + new_node_ratio*((u[end] - u[end-1])))
    insert!(TissueMech.k_init, idx, TissueMech.k_init[idx])
    insert!(TissueMech.k_max, idx, TissueMech.k_max[idx])
    insert!(TissueMech.λ, idx, TissueMech.λ[idx])
    insert!(TissueMech.t_0, idx, integrator.t)
    insert!(TissueMech.η, idx, TissueMech.η[idx])
    insert!(TissueMech.a, idx, TissueMech.a[idx])
    resize!(integrator, size(integrator.u, 1))
    nothing
end