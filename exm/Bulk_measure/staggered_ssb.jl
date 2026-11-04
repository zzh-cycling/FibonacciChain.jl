using FibonacciChain
using JLD2
using LinearAlgebra
using Random

"""
    simulate_staggered_ssb(L, λ_even, λ_odd, t; mode=:Born, seed=1,
                                        sign=true, samples=nothing)

Prepare |0> Π1|00…0> + |1> Πτ|00…0> with the original sector weights.
Each period applies odd sites 1:2:L with λ_odd, then even sites 2:2:L with λ_even.
Each event uses K_s = exp(-λ Π_i^s) / sqrt(1+exp(-2λ)) at its layer strength.
Here sign=true labels the suppressed vacuum channel, and false the τ channel.
Use mode=:Born to sample K_s, or :sample for a fixed record (default: all sign).
Time counts complete odd-then-even periods; no endpoint half-strength layer is added.
Fixed samples have size (2t,L÷2), with odd/even rows for odd/even bonds respectively.
Reference entropy is in nats and includes t=0.
"""
function simulate_staggered_ssb(
    L::Integer,
    λ_even::Real,
    λ_odd::Real,
    t::Integer;
    mode::Symbol = :Born,
    seed::Integer = 1,
    sign::Bool = true,
    samples::Union{Nothing,BitMatrix} = nothing,
)
    L >= 4 && iseven(L) || throw(ArgumentError("L must be even and at least 4"))
    all(λ -> isfinite(λ) && λ >= 0, (λ_even, λ_odd)) ||
        throw(ArgumentError("λ_even and λ_odd must be finite and nonnegative"))
    t >= 0 || throw(ArgumentError("t must be nonnegative"))
    mode ∈ (:Born, :sample) || throw(ArgumentError("mode must be :Born or :sample"))
    mode == :Born && samples !== nothing &&
        throw(ArgumentError("A fixed record requires mode=:sample"))
    if mode == :sample
        samples = isnothing(samples) ? BitMatrix(fill(sign, 2t, L÷2)) : copy(samples)
        size(samples) == (2t, L÷2) || throw(DimensionMismatch("samples must have size (2t,L÷2)"))
    end

    model = AnyonModel(FibonacciAnyon(), Int(L); pbc = true)
    basis = anyon_basis(model)
    initial_state = zeros(length(basis))
    initial_state[1] = 1
    ϕ = (1 + √5) / 2
    Ystate = topological_charge_operator(model) * initial_state
    state_y₁ = (Ystate .+ inv(ϕ) .* initial_state) ./ (ϕ + inv(ϕ))
    state_yτ = initial_state .- state_y₁
    joint_state = vcat(state_y₁, state_yτ)
    joint_state ./= norm(joint_state)

    rng = MersenneTwister(seed)
    configs = map((λ_odd, λ_even)) do λ
        MeasureConfig(
            τ = Float64(λ),
            t₂ = Int(t),
            mode = mode,
            rng = rng,
            enable_τ_eff = false,
        )
    end
    record = falses(2t, L÷2)
    sample_free_energy = zeros(2t)
    reference_entropy = zeros(t+1)
    sector_probabilities = zeros(t+1, 2)
    ρ = reference_rdm(model, [1], joint_state)
    reference_entropy[1] = ee(ρ)
    sector_probabilities[1, :] = real.(diag(ρ))

    for step in 1:t
        for layer in 1:2
            index = 2(step-1) + layer
            # The existing Kraus map favors its sign; complement it to suppress s.
            layer_sample = mode == :sample ? .!samples[index, :] : nothing
            # Library layer 2 is odd; layer 3 is even.
            outcome = reference_boundary_evolution(
                model, joint_state, configs[layer], layer_sample; layer_idx = index + 1,
            )
            joint_state = outcome.state
            all(isfinite, joint_state) || error("Invalid state at layer $index")
            record[index, :] = .!outcome.sample
            sample_free_energy[index] = outcome.free_energy
        end
        ρ = reference_rdm(model, [1], joint_state)
        reference_entropy[step+1] = ee(ρ)
        sector_probabilities[step+1, :] = real.(diag(ρ))
    end

    return (;
        L=Int(L), λ_even=Float64(λ_even), λ_odd=Float64(λ_odd),
        t=Int(t), times=collect(0:t), mode, seed=Int(seed), sign,
        measurement_sites_even=collect(2:2:L), measurement_sites_odd=collect(1:2:L),
        samples=record, sample_free_energy,
        initial_reference_entropy=first(reference_entropy), reference_entropy,
        sector_probabilities, reference_density_matrix=ρ, state=joint_state,
        initial_state="all_zero_fusion_path", reference_assignment="0:y=1; 1:y=tau",
        kraus_operator="exp(-lambda*Pi_i^s)/sqrt(1+exp(-2lambda))",
        sample_convention="true: suppressed vacuum channel; false: suppressed tau channel",
        time_unit="complete_odd_then_even_period", layer_order="odd_then_even", enable_τ_eff=false,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 5 || error(
        "Usage: julia --project exm/Bulk_measure/staggered_ssb.jl L lambda_even lambda_odd t Born|sample [seed] [sign] [output.jld2]",
    )
    L = parse(Int, ARGS[1])
    λ_even = parse(Float64, ARGS[2])
    λ_odd = parse(Float64, ARGS[3])
    t = parse(Int, ARGS[4])
    mode = Symbol(ARGS[5])
    seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 1
    sign = length(ARGS) >= 7 ? parse(Int, ARGS[7]) == 1 : true
    branch = mode == :Born ? "Born" : sign ? "sample_s1" : "sample_stau"
    output = length(ARGS) >= 8 ? ARGS[8] : joinpath(
        @__DIR__, "../data/Bulk_measure/staggered_ssb", branch,
        "L$L", "lambda_even$(λ_even)_lambda_odd$(λ_odd)", "t$t", "trajectory_seed$seed.jld2",
    )
    ispath(output) && error("Output already exists: $output")
    result = simulate_staggered_ssb(L, λ_even, λ_odd, t; mode, seed, sign)
    mkpath(dirname(abspath(output)))
    JLD2.jldsave(output; result...)
    println("Saved $output")
    println("S_R(0)=$(first(result.reference_entropy)), S_R(t)=$(last(result.reference_entropy))")
end
