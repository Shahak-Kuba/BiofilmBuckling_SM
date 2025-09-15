# User defined force function
function F(l, kₛ, a, restoring_force_type)
    if restoring_force_type == "hookean"
        return (kₛ .* (l .- a))'
    elseif restoring_force_type == "nonlinear"
        return (kₛ .* (1 ./ a .- 1 ./ l))'
    else
        return (kₛ .* a.^2 .* (1 ./ a .- 1 ./ l))'
    end
end


# restoring force equation
#Fs = (k, a, u1, u2) -> k .* (abs.(u1 .- u2) .- a)

# spring stiffness evolution equation (logistic growth)
#k_log_growth = (K_init, K_max, λ, t_0, t_curr) -> (K_init.*K_max.*exp.(λ.*(t_curr .- t_0))) ./ ((K_max .- K_init) .+ K_init .* exp.(λ .* (t_curr .- t_0)))

"""
    δ(rᵢ₊₁, rᵢ)

Calculate the Euclidean distance between two points `rᵢ₊₁` and `rᵢ`.

# Arguments
- `rᵢ₊₁`: The first point in space.
- `rᵢ`: The second point in space.

# Returns
The Euclidean distance between the two points.
"""
#δ(rᵢ₊₁, rᵢ) = .√(sum((rᵢ₊₁ - rᵢ).^2,dims=size(rᵢ)))

function δ(rᵢ₊₁, rᵢ)
    val = zeros(1,size(rᵢ,2))
    val .= .√(sum((rᵢ₊₁- rᵢ).^2,dims=1))
    return val
end

"""
    τ(rᵢ₊₁, rᵢ₋₁)

Calculate the unit tangent vector between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.

# Returns
The unit tangent vector between the two points.
"""
τ(rᵢ₊₁, rᵢ₋₁) = (rᵢ₊₁ - rᵢ₋₁) ./ δ(rᵢ₊₁, rᵢ₋₁)

"""
    n(rᵢ₊₁, rᵢ₋₁, type)

Calculate the unit normal vector at a point `rᵢ` between two neighboring points `rᵢ₊₁` and `rᵢ₋₁`. The orientation of the normal vector depends on the specified `type`.

# Arguments
- `rᵢ₊₁`: The point after the central point in space.
- `rᵢ₋₁`: The point before the central point in space.
- `type`: A string specifying the orientation of the normal vector, either "inward" or any other value for outward orientation.

# Returns
The unit normal vector at the point.
"""
function n(rᵢ₊₁, rᵢ₋₁,type) 
    if type == "inward"
        -oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    else
        oftype(τ(rᵢ₊₁, rᵢ₋₁),vcat(transpose.([-τ(rᵢ₊₁, rᵢ₋₁)[:,2], τ(rᵢ₊₁, rᵢ₋₁)[:,1]])...)')
    end
end


function calc_l_ρ_dv_τ_n(rᵢ₊₁, rᵢ, rᵢ₋₁, growth_dir)
    l = .√(sum((rᵢ₊₁ - rᵢ).^2, dims=1)) # calculating length
    ρ = 1.0 ./ l # calculating density
    dv = (rᵢ₊₁ - rᵢ) ./ l # directional vector
    τ = (rᵢ₊₁ - rᵢ₋₁) ./ .√(sum((rᵢ₊₁ - rᵢ₋₁).^2, dims=1)) # tangent vector

    if growth_dir == "inward" # normal vector
        n = [-1.0, 1.0] .* circshift(dv, 1)
    else
        n = [1.0, -1.0] .* circshift(dv, 1)
    end

    return l, ρ, dv, τ, n # returning all calculated values
end

