include("config.jl")
include("../reference_correlations.jl")

"""Two-reference MI in the exact OBF circuit, starting from the even-parity |+> state."""
function obf_reference_sample(L::Int, λ::Float64, τ_idx::Int, separation::Int, seed::Int;
    warmup::Int, final_periods::Int, record = nothing,
)
    model = obf_model(L, λ)
    state = normalize(ones(ComplexF64, length(anyon_basis(model))))
    result = reference_correlation_sample(model, state; τ = τlis[τ_idx], seed,
        warmup, separation, final_periods, record,
        x1 = separation == 0 ? 1 : L ÷ 2 + 1, x2 = L ÷ 2 + 1)
    return (; result..., λ, τ_idx, initial_state = "plus_product", circuit = "OBF")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 9 || error("Usage: julia --project=. exm/Born_OBF/corr_calculate.jl L lambda tau_idx warmup final_periods dt_csv seed_start seed_end output_dir")
    L = parse(Int, ARGS[1])
    λ = parse(Float64, ARGS[2])
    τ_idx = parse(Int, ARGS[3])
    warmup, final_periods = parse.(Int, ARGS[4:5])
    separations = parse.(Int, split(ARGS[6], ','))
    seeds = parse(Int, ARGS[7]):parse(Int, ARGS[8])
    generate = (dt, seed) -> obf_reference_sample(L, λ, τ_idx, dt, seed; warmup, final_periods)
    result = run_reference_correlation_scan(generate, L, separations, seeds, ARGS[9])
    println("alpha = $(result.alpha), bootstrap SEM = $(result.alpha_sem)")
end
