import BiofilmBucklingSM as BBSM
using DifferentialEquations
using CairoMakie
using ElasticArrays

N = 30; m = 1;
M = N * m; # total number of nodes

T1 = BBSM.TissueMechProperties_t(k=5.0, a=3.0, η=5, restoring_force="Hookean")
TissueMech = BBSM.generate_homogeneous_population(T1, N, m)
SimTime = BBSM.SimTime_t(Tmax=20.0, δt=0.01, event_δt=2.0)

# solving problem
#event_cb = PeriodicCallback(TissueRelaxation1D.tissue_dep_affect!, SimTime.event_δt; save_positions = (false, false))
#cbs = CallbackSet(event_cb)

# setup initial condition
u0 = ElasticMatrix{Float64}(undef,2,M)
u0[1,:] = collect(LinRange(0,75, M))
u0[2,:] .= 0.0
#u0[2,7] = 0.025 # perturb the middle node
u0[2,15] = 0.05 # perturb the middle node
#u0[2,23] = 0.025 # perturb the middle node


prob, p = BBSM.SetupODEproblem(TissueMech, SimTime, u0);

sol = solve(prob, Euler(), dt = SimTime.δt, dtmax = SimTime.δt)#, callback=cbs);

f = Figure(size=(650,500));
ax = Axis(f[1,1], aspect=1, xlabel=L"$x\; \text{[μm]}$", ylabel=L"$\text{Time [days]}$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24);
for ii in 1:500:size(sol.u,1)
    x = sol.u[ii][1,:]
    y = sol.u[ii][2,:]
    lines!(ax, x, y, color=sol.t[ii], colormap=:viridis, linewidth=3, colorrange=(0,SimTime.Tmax));
end
Colorbar(f[1,2], label = L"$\text{Time [days]}$", labelsize=32, ticklabelsize=24, colorrange=(0,SimTime.Tmax), colormap = :viridis);

display(f)

using CairoMakie

frames = [Matrix(sol.u[ii]) for ii in 1:size(sol.u,1)]

f = Figure(size=(650,500))
ax = Axis(f[1,1], aspect=1, xlabel=L"$x\; \text{[μm]}$", ylabel=L"$y$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24)
lineplot = lines!(ax, frames[1][1,:], frames[1][2,:], color=:blue, linewidth=3)
display(f)

record(f, "biofilm_evolution.mp4", 1:length(frames); framerate=100) do i
    lineplot.arg1.value.x = frames[i][1,:]
    lineplot.arg2.value.x = frames[i][2,:]
    ax.title = "Time = $(round(sol.t[i], digits=2))"
end


"""
    animate_xy_frames(frames; filename="motion.mp4", framerate=30, markersize=8, title=nothing)

Create an animation from `frames::Vector{<:AbstractMatrix}` where each element is a 2×N matrix:
- row 1 = x coordinates
- row 2 = y coordinates

Saves an .mp4 (or .gif if you pass a .gif filename) and returns the Figure.
"""
function animate_xy_frames(frames::AbstractVector{<:AbstractMatrix};
                           filename::AbstractString = "motion.mp4",
                           framerate::Integer = 30,
                           markersize = 8,
                           title = nothing)

    @assert !isempty(frames) "frames is empty"
    # Validate shapes and collect global bounds to keep axes stable
    xmin = Inf; xmax = -Inf; ymin = Inf; ymax = -Inf
    npts = nothing
    for (k, M) in enumerate(frames)
        @assert size(M, 1) == 2 "Frame $k must be 2×N (got size $(size(M)))"
        npts === nothing && (npts = size(M, 2))
        @assert size(M, 2) == npts "All frames must have the same number of points (frame $k differs)"
        xk, yk = @view(M[1, :]), @view(M[2, :])
        xmin = min(xmin, minimum(xk)); xmax = max(xmax, maximum(xk))
        ymin = min(ymin, minimum(yk)); ymax = max(ymax, maximum(yk))
    end
    # Pad bounds a touch so points near the edge aren't clipped
    dx = max(1e-9, 0.05 * max(1e-9, xmax - xmin))
    dy = max(1e-9, 0.05 * max(1e-9, ymax - ymin))
    xlims = (xmin - dx, xmax + dx)
    ylims = (ymin - dy, ymax + dy)

    CairoMakie.activate!()  # ensure we're using CairoMakie backend

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], title = something(title, "XY animation"),
              xlabel = "x", ylabel = "y")

    # Observables for positions
    xs = Observable(frames[1][1, :])
    ys = Observable(frames[1][2, :])

    # Plot: scatter; change to lines! if your columns represent a polyline
    scatter!(ax, xs, ys; markersize)

    # Optional: if you want connecting lines between points in each frame, uncomment:
    # lines!(ax, xs, ys; transparency = true)

    display(fig)  # helpful in notebooks

    record(fig, filename, 1:length(frames); framerate = framerate) do i
        xs[] = frames[i][1, :]
        ys[] = frames[i][2, :]
    end

    println("Saved animation to: $filename")
    return fig
end

animate_xy_frames(frames; filename="swarm.mp4", framerate=30, markersize=6)
