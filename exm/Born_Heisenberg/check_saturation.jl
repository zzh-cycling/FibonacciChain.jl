# Check saturation of the half-chain EE and free-energy dynamics for a given
# τ index, using the collected ensemble files. Prints the plateau value and
# the saturation time (in units of L) for each Δ, defined as the last time
# the trajectory deviates from the plateau by more than `tol`.
#
# Usage: julia --project check_saturation.jl [τ_idx] [tol]
using JLD2
using Statistics
using Printf

include(joinpath(@__DIR__, "config.jl"))

ind = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10
tol = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.02
τ = τlis[ind]

function satinfo(x::AbstractVector, tol::Float64; nskip_end::Int = 2)
    # Exclude the last nskip_end periods: they use the τ/2 "effective" layer
    # (enable_τ_eff) and are not part of the steady-state trajectory.
    n = length(x) - nskip_end
    xx = @view x[1:n]
    plateau = mean(xx[end-div(n, 8)+1:end])
    thr = max(tol * abs(plateau), 1e-3)
    bad = findall(v -> abs(v - plateau) > thr, xx)
    return plateau, isempty(bad) ? 1 : maximum(bad) + 1
end

for L in [8, 10, 12, 14]
    println("==== L = $L, τ = $τ (ind = $ind), tol = $(100*tol)% (last 2 periods excluded) ====")
    for Δ in Δlis
        t, _, _ = get_born_dynamics_params(ind, L, Δ)
        f = "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Observables_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2"
        if !isfile(f)
            @printf("Δ=%5.1f  MISSING %s\n", Δ, f)
            continue
        end
        d = load(f)
        ee = d["average_EE_tlis"]
        fe = sum(reshape(d["time_FElis"], 2, :), dims = 1)[:]
        ee_plateau, t_ee = satinfo(ee, tol)
        fe_plateau, t_fe = satinfo(fe, tol)
        _, t_ee5 = satinfo(ee, 0.05)
        _, t_fe5 = satinfo(fe, 0.05)
        @printf(
            "Δ=%5.1f  ee_plateau=%8.4f  t_sat(ee)/L=%6.2f (5%%: %5.2f)  fe_plateau=%9.5f  t_sat(fe)/L=%6.2f (5%%: %5.2f)\n",
            Δ, ee_plateau, t_ee / L, t_ee5 / L, fe_plateau, t_fe / L, t_fe5 / L
        )
    end
end

# Detailed ee curves at integer multiples of L for a few representative Δ
println("\n==== EE(t/L) curves ====")
for L in [8, 12, 14]
    for Δ in [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0]
        t, _, _ = get_born_dynamics_params(ind, L, Δ)
        f = "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Observables_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2"
        isfile(f) || continue
        ee = load(f, "average_EE_tlis")
        vals = [@sprintf("%.4f", ee[k*L]) for k in 1:t]
        println("L=$L Δ=$Δ: " * join(vals, " "))
    end
end

# Head/tail of the trajectories: check for genuine transients vs the τ/2
# artifact in the last layer (enable_τ_eff).
println("\n==== head/tail of ee and fe (per-period) ====")
for (L, Δ) in [(14, 0.7), (14, 0.9), (10, 0.7), (8, 0.7), (12, 0.9), (8, 0.6)]
    t, _, _ = get_born_dynamics_params(ind, L, Δ)
    f = "exm/data/Heisenberg/Born_dynamics_records/L$(L)/τ$(τ)/Observables_L$(L)_τ$(τ)_Δ$(Δ)_t$(t).jld2"
    isfile(f) || continue
    d = load(f)
    ee = d["average_EE_tlis"]
    fe = sum(reshape(d["time_FElis"], 2, :), dims = 1)[:]
    n = length(ee)
    fmt(v) = join([@sprintf("%.4f", x) for x in v], " ")
    println("L=$L Δ=$Δ n=$n")
    println("  ee[1:$(2L)]   = ", fmt(ee[1:2L]))
    println("  ee[end-5:end] = ", fmt(ee[end-5:end]))
    println("  fe[1:$(2L)]   = ", fmt(fe[1:2L]))
    println("  fe[end-5:end] = ", fmt(fe[end-5:end]))
end
