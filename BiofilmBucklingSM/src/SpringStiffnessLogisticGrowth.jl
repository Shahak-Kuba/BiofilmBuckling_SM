
k_log_growth = (K_init, K_max, λ, t) -> (K_init.*K_max.*exp.(λ.*t)) ./ ((K_max .- K_init) .+ K_init .* exp.(λ .* t))
k_log_growth = (K_init, K_max, λ, t_0, t_curr) -> (K_init.*K_max.*exp.(λ.*(t_curr .- t_0))) ./ ((K_max .- K_init) .+ K_init .* exp.(λ .* (t_curr .- t_0)))

k_max = 100;
k_init = 5;
λ = 1.0;

t = 0:0.1:10;
k_values_0 = k_log_growth(k_init, k_max, λ, 0, t);
using CairoMakie
f = Figure(size = (350,350))
ax = Axis(f[1,1], aspect=1, xlabel=L"$\text{Time [days]}$", ylabel=L"$k(t)\; \text{[μN/μm]}$")
lines!(ax, t, k_values_0, color=:blue, linewidth=3)
display(f)
save("SpringStiffnessLogisticGrowth.png", f, px_per_unit=2.0)