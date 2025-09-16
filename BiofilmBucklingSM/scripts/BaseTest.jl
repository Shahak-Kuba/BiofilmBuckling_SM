import BiofilmBucklingSM as BBSM
using DifferentialEquations
using CairoMakie
using ElasticArrays

N = 65; m = 1;
M = N * m; # total number of nodes

T1 = BBSM.TissueMechProperties_t(k=10.0, a=1, η=5, restoring_force="Hookean")
TissueMech = BBSM.generate_homogeneous_population(T1, N, m)
SimTime = BBSM.SimTime_t(Tmax=40.0, δt=0.01, event_δt=2.0)

# solving problem
#event_cb = PeriodicCallback(TissueRelaxation1D.tissue_dep_affect!, SimTime.event_δt; save_positions = (false, false))
#cbs = CallbackSet(event_cb)

# substrate function
S = (x) -> 0.05 .* sin.((π/20) .* x) 

# setup initial condition
Bell_curve = (x) -> 0.5 .* exp.(-((x .- 30).^2) ./ 20)
u0 = ElasticMatrix{Float64}(undef,2,M)
u0[1,:] = collect(LinRange(0,60, M))
u0[2,:] = Bell_curve.(u0[1,:]) .+ 0.06


prob, p = BBSM.SetupODEproblem(TissueMech, SimTime, u0, S);

sol = solve(prob, Euler(), dt = SimTime.δt, dtmax = SimTime.δt)#, callback=cbs);

f = Figure(size=(650,500));
ax = Axis(f[1,1], aspect=1, xlabel=L"$x\; \text{[μm]}$", ylabel=L"$\text{Time [days]}$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24);
lines!(ax, collect(LinRange(0,60, 200)), S.(collect(LinRange(0,60, 200))), color=:black, linewidth=5, label="Substrate");
for ii in 1:500:size(sol.u,1)
    x = sol.u[ii][1,:]
    y = sol.u[ii][2,:]
    lines!(ax, x, y, color=sol.t[ii], colormap=:viridis, linewidth=3, colorrange=(0,SimTime.Tmax));
end
Colorbar(f[1,2], label = L"$\text{Time [days]}$", labelsize=32, ticklabelsize=24, colorrange=(0,SimTime.Tmax), colormap = :viridis);

display(f)

# Animation

frames = [Matrix(sol.u[ii]) for ii in 1:size(sol.u,1)]

f = Figure(size=(650,500))
ax = Axis(f[1,1], aspect=1, xlabel=L"$x\; \text{[μm]}$", ylabel=L"$y$", xlabelsize=32, ylabelsize=32, xticklabelsize=24, yticklabelsize=24, limits=((-0.1,60),(-0.1,5)));
substrate_line = lines!(ax, collect(LinRange(0,60, 200)), S.(collect(LinRange(0,60, 200))), color=:black, linewidth=5, label="Substrate");
lineplot = lines!(ax, frames[1][1,:], frames[1][2,:], color=:blue, linewidth=3)
scatterplot = scatter!(ax, frames[1][1,:], frames[1][2,:], color=:black, markersize=8)
display(f)

record(f, "biofilm_evolution.mp4", 1:length(frames); framerate=60) do i
    lineplot[1] = frames[i][1,:]
    lineplot[2] = frames[i][2,:]
    scatterplot[1] = frames[i][1,:]
    scatterplot[2] = frames[i][2,:]
    ax.title = "Time = $(round(sol.t[i], digits=2))"
end

