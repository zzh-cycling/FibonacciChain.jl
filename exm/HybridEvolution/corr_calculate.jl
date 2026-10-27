include("../reference_correlations.jl")

"""Two-reference MI in the exclusive hybrid circuit, starting from the all-zero fusion path."""
function hybrid_reference_sample(L::Int, p::Float64, separation::Int, seed::Int;
    warmup::Int, final_periods::Int, τ::Float64 = Inf, record = nothing,
)
    model = AnyonModel(FibonacciAnyon(), L; pbc = true)
    state = zeros(ComplexF64, length(anyon_basis(model)))
    state[1] = 1
    config = HybridConfig(; τ, p, t₂ = warmup + separation + final_periods,
        random_angles = true, enable_τ_eff = false)
    result = reference_correlation_sample(model, state; τ, seed,
        warmup, separation, final_periods, hybrid = config, record,
        x1 = separation == 0 ? 1 : L ÷ 2 + 1, x2 = L ÷ 2 + 1)
    return (; result..., p, random_angles = true,
        initial_state = "all_zero_fusion_path", circuit = "HybridEvolution")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 8 || error("Usage: julia --project=. exm/HybridEvolution/corr_calculate.jl L p warmup final_periods dt_csv seed_start seed_end output_dir")
    L = parse(Int, ARGS[1])
    p = parse(Float64, ARGS[2])
    warmup, final_periods = parse.(Int, ARGS[3:4])
    separations = parse.(Int, split(ARGS[5], ','))
    seeds = parse(Int, ARGS[6]):parse(Int, ARGS[7])
    generate = (dt, seed) -> hybrid_reference_sample(L, p, dt, seed; warmup, final_periods)
    result = run_reference_correlation_scan(generate, L, separations, seeds, ARGS[8])
    println("alpha = $(result.alpha), bootstrap SEM = $(result.alpha_sem)")
end
