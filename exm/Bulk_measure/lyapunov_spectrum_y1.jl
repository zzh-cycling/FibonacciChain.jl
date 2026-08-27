using FibonacciChain
using JLD2
using LinearAlgebra
using Random

"""
    initial_y1_state(model) -> Vector

Project the all-`τ` fusion path into the `y=1` sector and normalize it. This is
the state used to sample the Born trajectory. In the repository normalization
the `y=1` sector has topological charge eigenvalue `ϕ`, while the `y=τ` sector
has eigenvalue `-1/ϕ`, so `P1 = (Y + ϕ⁻¹ I) / (ϕ + ϕ⁻¹)`.
"""
function initial_y1_state(model::AnyonModel{FibonacciAnyon})
    model.pbc || error("The topological y=1 sector requires periodic boundaries")
    ϕ = (1 + √5) / 2
    Y = topological_charge_operator(model)
    P1 = (Y + inv(ϕ) * I(size(Y, 1))) / (ϕ + inv(ϕ))
    state = zeros(Float64, length(anyon_basis(model)))
    state[1] = 1.0
    projected = P1 * state
    weight = real(dot(projected, projected))
    weight > 1e-14 || error("The all-τ path has zero y=1 weight")
    return projected / sqrt(weight)
end

"""
    simulate_y1_lyapunov(L, τ, periods; n_states=10, trajectory_seed=1)

Generate a Born trajectory from `P1|τ⋯τ⟩` (exact state-vector evolution; the
topological Y operator has no MPS representation, so no MPS variant exists),
then compute the Lyapunov spectrum of the corresponding transfer-matrix
product restricted to the full `y=1` sector via `lyapunov_spectrum_topological_sector`.
"""
function simulate_y1_lyapunov(
    L::Int,
    τ::Real,
    periods::Int;
    n_states::Int = 10,
    trajectory_seed::Int = 1,
)
    periods >= 1 || throw(ArgumentError("periods must be positive"))
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    initial_state = initial_y1_state(model)
    born_config = MeasureConfig(
        τ = Float64(τ),
        t₂ = periods,
        mode = :Born,
        rng = MersenneTwister(trajectory_seed),
        enable_τ_eff = false,
        track_y_expectation = true,
    )
    trajectory = bulk_evolution(model, initial_state, born_config)
    spectrum = lyapunov_spectrum_topological_sector(
        model,
        Float64(τ),
        trajectory.samples;
        sector = :trivial,
        n_states = n_states,
    )

    ϕ = (1 + √5) / 2
    maximum(abs.(trajectory.y_expectation_values .- ϕ)) < 1e-5 ||
        error("the Born trajectory left the y=1 sector")
    maximum(spectrum.sector_leakage) < 1e-9 ||
        error("the Lyapunov frame left the y=1 sector")
    return trajectory, spectrum
end


function main(args = ARGS)
    length(args) in (5, 6) || error(
        "Usage: julia --project=. exm/Bulk_measure/lyapunov_spectrum_y1.jl " *
        "L tau periods n_states seed [output.jld2]",
    )
    L = parse(Int, args[1])
    τ = parse(Float64, args[2])
    periods = parse(Int, args[3])
    n_states = parse(Int, args[4])
    seed = parse(Int, args[5])

    trajectory, spectrum = simulate_y1_lyapunov(
        L,
        τ,
        periods;
        n_states = n_states,
        trajectory_seed = seed,
    )
    output_path = if length(args) == 6
        args[6]
    else
        joinpath(
            @__DIR__,
            "..",
            "data",
            "Bulk_measure",
            "lyapunov_spectrum_y1",
            "L$(L)_tau$(τ)_periods$(periods)_seed$(seed).jld2",
        )
    end
    mkpath(dirname(output_path))
    jldsave(
        output_path;
        L = L,
        τ = τ,
        periods = periods,
        n_states = n_states,
        trajectory_seed = seed,
        topological_sector = "y=1",
        y_eigenvalue = (1 + √5) / 2,
        sector_dimension = spectrum.sector_dimension,
        samples = trajectory.samples,
        local_log_stretches = spectrum.local_log_stretches,
        lyapunov_exponents = spectrum.lyapunov_exponents,
        free_energy_spectrum = spectrum.free_energy_spectrum,
        sector_leakage = spectrum.sector_leakage,
        y_expectation_values = trajectory.y_expectation_values,
    )

    println("saved: $output_path")
    println("y=1 sector dimension: $(spectrum.sector_dimension)")
    println("max sector leakage: $(maximum(spectrum.sector_leakage))")
    println("final Lyapunov exponents: $(spectrum.lyapunov_exponents[:, end])")
    println("final free-energy spectrum: $(spectrum.free_energy_spectrum[:, end])")
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
