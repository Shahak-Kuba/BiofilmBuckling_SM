
# user specified substrate function to calculate y coords of substrate surface
y = (x) -> x^3

function substrate_surface(y,x)
    return [x; y(x)]
end

function substrate_tangents(y,x,type)
    dx = 1e-5
    dy = y(x + dx) - y(x)
    τ = [dx; dy] ./ √(dx^2 + dy^2)
    return τ
end

# calculate normal vectors of substrate surface
function substrate_normal(y,x,type)
    dx = 1e-5
    dy = y(x + dx) - y(x)
    τ = [dx; dy] ./ √(dx^2 + dy^2)
    if type == "inward"
        n = [-τ[2]; τ[1]]
    else
        n = [τ[2]; -τ[1]]
    end
    return n
end

# function to project position derivative vector onto substrate normal direction
function proj_dr_onto_S_normal(dr, y, x, type)
    n = substrate_normal(y, x, type)
    return (dr' * n)[1] * n
end

function check_dr(r, dr, y, dt)
    # checking that derivative does not move node below substrate
    r_next = r .+ (dr .* dt)
    x_next = r_next[1,:]
    for ii in eachindex(x_next)
        #println("x_next: ", x_next[ii], " y(x_next): ", y(x_next[ii]), " r_next[2,ii]: ", r_next[2,ii])
        if r_next[2,ii] < y(x_next[ii])
            dr[:,ii] = dr[:,ii] - proj_dr_onto_S_normal(dr[:,ii], y, x_next[ii], "inward")
        end
    end

    # checking that derivatives do no move nodes to the left or right of neighbors and ensureing there is a minimum distance between nodes ii-2 and ii, and ii and ii+2
    
    """
    M = size(r,2)
    for ii in [2,M-1]
        if norm(r_next[1,ii] - r_next[1,ii-1]) < 2 
            if dr[1,ii] < 0
                dr[1,ii] = 0.1 * dr[1,ii]
            end
        end
        if norm(r_next[1,ii] - r_next[1,ii+1]) < 2
            if dr[1,ii] > 0
                dr[1,ii] = 0.1 * dr[1,ii]
            end
        end
    end
    for ii in 3:M-2
        if norm(r_next[1,ii] - r_next[1,ii-2]) < 2 
            if dr[1,ii] < 0
                dr[1,ii] = 0.1 * dr[1,ii]
            end
        end
        if norm(r_next[1,ii] - r_next[1,ii+2]) < 2
            if dr[1,ii] > 0
                dr[1,ii] = 0.1 * dr[1,ii]
            end
        end
    end
    """

    return dr
end





















 using CairoMakie

# Define x range and sample points for normals
x_vals = range(-2, 2, length=200)
x_normals = range(-2, 2, length=10)

# Compute substrate surface points
surface_points = [substrate_surface(y, x) for x in x_vals]
xs = [p[1] for p in surface_points]
ys = [p[2] for p in surface_points]

# Compute normal vectors at selected points
normal_points = [substrate_surface(y, x) for x in x_normals]
tangents = [substrate_tangents(y, x, "inward") for x in x_normals]
normals = [substrate_normal(y, x, "inward") for x in x_normals]

fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], xlabel="x", ylabel="y", title="Substrate Surface & Normals")

lines!(ax, xs, ys, color=:blue, linewidth=2)

normal_xs = [p[1] for p in normal_points]
normal_ys = [p[2] for p in normal_points]
tangents_us = [τ[1] for τ in tangents]
tangents_vs = [τ[2] for τ in tangents]
normals_us = [n[1] for n in normals]
normals_vs = [n[2] for n in normals]
arrows2d!(ax, normal_xs, normal_ys, tangents_us, tangents_vs, color=:green)
arrows2d!(ax, normal_xs, normal_ys, normals_us, normals_vs, color=:red)

fig