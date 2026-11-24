using FibonacciChain, ITensors, ITensorMPS, LinearAlgebra, Random

function _staggered_mps_branches(joint, ref)
    [MPS(vcat([joint[1] * onehot(ref => r) * joint[2]], joint[3:end])) for r in 1:2]
end

function _staggered_mps_sector_state(ψ, Y, reference; cutoff, maxdim)
    ϕ = (1 + sqrt(5)) / 2
    yψ = apply(Y, ψ; cutoff, maxdim)
    reference == 1 && return add(yψ, inv(ϕ) * ψ; cutoff, maxdim) / sqrt(5)
    add(ϕ * ψ, -yψ; cutoff, maxdim) / sqrt(5)
end

function _staggered_mps_region_entropy(ψ, region)
    grouped = movesites(ψ, region, collect(1:length(region)); cutoff=1e-14)
    ee_mps(grouped, length(region))
end

function _staggered_mps_final(joint, sites, Y, ref; cutoff, maxdim)
    L = length(sites)
    width = L ÷ 4
    A, B, C = collect(1:width), collect(width+1:2width), collect(2width+1:3width)
    regions = [A, B, C, vcat(A,B), vcat(A,C), vcat(B,C), vcat(A,B,C)]
    branches = _staggered_mps_branches(joint, ref)
    weights = [real(inner(ψ, ψ)) for ψ in branches]
    y = fill(NaN, 2)
    profiles = fill(NaN, L-1, 2)
    entropies = fill(NaN, 7, 2)
    correlation = zeros(L, L)
    occupation_weight = 0.0
    for (r, ψ) in enumerate(branches)
        weights[r] > 0 || continue
        rmul!(ψ, inv(sqrt(weights[r])))
        y[r] = real(inner(ψ', Y, ψ))
        profiles[:,r] = [ee_mps(ψ, x) for x in 1:L-1]
        entropies[:,r] = [_staggered_mps_region_entropy(ψ, region) for region in regions]
        τψ = _staggered_mps_sector_state(ψ, Y, 2; cutoff, maxdim)
        occupation_weight += weights[r] * real(inner(τψ, τψ))
        frame = map(sites) do site
            occupied = apply(op("Proj1", site), τψ; cutoff=0.0)
            _staggered_mps_sector_state(occupied, Y, 1; cutoff, maxdim)
        end
        for i in 1:L, j in 1:i
            correlation[i,j] += weights[r] * real(inner(frame[i], frame[j]))
        end
    end
    correlation = Matrix(Hermitian(correlation, :L))
    I3 = vec(sum(entropies[1:3,:]; dims=1) - sum(entropies[4:6,:]; dims=1) + entropies[7:7,:])
    (; final_subsystem_entropy=profiles, I3, I3_component_entropies=entropies,
        y_expectation=sum(weights[r] * y[r] for r in 1:2 if weights[r] > 0),
        y_expectation_by_reference=y, reference_weights=weights,
        occupation_correlation_tau=correlation, occupation_weight_tau=occupation_weight)
end

"""
    simulate_staggered_ssb_mps(L, λ_even, λ_odd, t; seed=1, outcomes=nothing,
                               cutoff=1e-12, maxdim=128)

Born sample or replay the periodic staggered circuit, applying the odd layer
first. The reference basis is |0⟩ for y=1 and |1⟩ for y=τ. Outcome bits label
the suppressed vacuum (true) or τ (false) channel, as in the exact code.
Born probabilities are evaluated before truncating the selected branch.
All observables are returned in memory. Occupation projections use twice the
evolution bond-dimension cap. Compare bond dimensions and cutoffs before treating the trajectory as a converged large-system calculation.
"""
function simulate_staggered_ssb_mps(
    L::Int,
    λ_even::Real,
    λ_odd::Real,
    t::Int;
    seed::Int=1,
    outcomes::Union{Nothing,BitMatrix}=nothing,
    cutoff::Float64=1e-12,
    maxdim::Int=128,
)
    iseven(L) && L >= 4 || throw(ArgumentError("L must be even and at least 4"))
    t >= 0 || throw(ArgumentError("t must be nonnegative"))
    all(λ -> !isnan(λ) && λ >= 0, (λ_even, λ_odd)) ||
        throw(ArgumentError("layer strengths must be nonnegative"))
    maxdim >= 1 && cutoff >= 0 || throw(ArgumentError("invalid truncation settings"))
    isnothing(outcomes) || size(outcomes) == (2t, L÷2) ||
        throw(DimensionMismatch("outcomes must have size (2t,L÷2)"))
    rng = MersenneTwister(seed)
    model = AnyonModel(FibonacciAnyon(), L; pbc=true)
    ψ, sites = initial_mps(L)
    Y = topological_charge_mpo(sites; pbc=true)
    ref = Index(2, "Qubit,Site,n=ref")
    joint_sites = vcat([ref], sites)
    y1, yτ = (_staggered_mps_sector_state(ψ, Y, r; cutoff, maxdim) for r in 1:2)
    joint = add(FibonacciChain._embed_system_state_with_ancilla(y1, ref, 1),
        FibonacciChain._embed_system_state_with_ancilla(yτ, ref, 2); cutoff, maxdim)
    rmul!(joint, inv(sqrt(real(inner(joint, joint)))))

    operators = map(enumerate((λ_odd, λ_even))) do (layer, λ)
        sites_to_measure, layer_model, strength = FibonacciChain._obtain_measurement_config(
            model, layer + 1, min(Float64(λ), 1000.0))
        [(site, FibonacciChain._joint_measurement_operator_application(
            layer_model, sites, joint_sites, site, strength, false),
          FibonacciChain._joint_measurement_operator_application(
            layer_model, sites, joint_sites, site, strength, true)) for site in sites_to_measure]
    end
    record = falses(2t, L÷2)
    reference_entropy = zeros(t+1)
    sector_probabilities = zeros(t+1, 2)
    S_half = fill(NaN, t+1, 2)
    sample_free_energy = zeros(2t)
    max_bond_dimension = maxlinkdim(joint)
    max_discarded_weight = 0.0
    ρ = reference_rdm(model, joint)
    for period in 0:t
        if period > 0
            for layer in 1:2
                index = 2(period-1) + layer
                for (event, (site, op_false, op_true)) in enumerate(operators[layer])
                    if isnothing(outcomes)
                        candidate = apply(op_false, joint; cutoff=0.0)
                        p_false = clamp(real(inner(candidate, candidate)), 0.0, 1.0)
                        suppressed_vacuum = rand(rng) < p_false
                        if !suppressed_vacuum
                            candidate = apply(op_true, joint; cutoff=0.0)
                        end
                    else
                        suppressed_vacuum = outcomes[index,event]
                        candidate = apply(suppressed_vacuum ? op_false : op_true, joint; cutoff=0.0)
                    end
                    probability = real(inner(candidate, candidate))
                    isfinite(probability) && probability > 0 ||
                        error("invalid branch probability at layer $index, event $event")
                    sample_free_energy[index] -= log(probability)
                    record[index,event] = suppressed_vacuum
                    rmul!(candidate, inv(sqrt(probability)))
                    # Only the two bonds inside a bulk gate can grow. The
                    # wrap-around MPO can grow bonds across the whole chain.
                    support = 1 < site < L ? (site:site+2) : (1:L+1)
                    truncate!(candidate; cutoff, maxdim, site_range=support)
                    retained_weight = real(inner(candidate, candidate))
                    discarded_weight = max(0.0, 1 - retained_weight)
                    max_discarded_weight = max(max_discarded_weight, discarded_weight)
                    rmul!(candidate, inv(sqrt(retained_weight)))
                    joint = candidate
                    max_bond_dimension = max(max_bond_dimension, maxlinkdim(joint))
                end
            end
            ρ = reference_rdm(model, joint)
        end
        ρ ./= real(tr(ρ))
        reference_entropy[period+1] = ee(ρ)
        sector_probabilities[period+1,:] = real.(diag(ρ))
        for (r, branch) in enumerate(_staggered_mps_branches(joint, ref))
            weight = real(inner(branch, branch))
            weight > 0 || continue
            rmul!(branch, inv(sqrt(weight)))
            S_half[period+1,r] = ee_mps(branch, L÷2)
        end
    end

    free_energy = vcat(0.0, cumsum(sample_free_energy[1:2:end] + sample_free_energy[2:2:end]))
    (; seed, reference_entropy, sector_probabilities, S_half, sample_free_energy, free_energy,
        _staggered_mps_final(joint, sites, Y, ref; cutoff, maxdim=2maxdim)..., reference_density_matrix=ρ,
        outcomes=record, max_bond_dimension, max_discarded_weight)
end
