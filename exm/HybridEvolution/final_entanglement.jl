"""Four consecutive, nearly equal blocks; equal quarters when L is divisible by four."""
function final_entanglement_partition(L::Int)
    L >= 4 || throw(ArgumentError("Four nonempty blocks require L >= 4"))
    edges = [fld(k * L, 4) for k in 0:4]
    return [collect((edges[k] + 1):edges[k + 1]) for k in 1:4]
end

"""Final prefix EE profile and I3 in the fusion-path site convention, using natural logs.

`I3_entropies` is ordered as A, B, C, AB, AC, BC, ABC. The fourth block D
completes the chain. The input is a pure state, so complements can be used
to keep exact reduced density matrices small.
"""
function final_entanglement(model::AnyonModel, state::Vector)
    L = model.N
    I3_partition = final_entanglement_partition(L)
    A, B, C, D = I3_partition
    entropy(sites) = ee(anyon_rdm(model,
        length(sites) <= L ÷ 2 ? sites : setdiff(collect(1:L), sites), state))
    subsystem_sizes = collect(1:(L - 1))
    S_subsystem_final = anyon_eelis(model, state)
    I3_entropies = [S_subsystem_final[last(A)], entropy(B), entropy(C),
        S_subsystem_final[last(B)], entropy(vcat(A, C)), entropy(vcat(B, C)),
        S_subsystem_final[last(C)]]
    I3_final = sum(I3_entropies[1:3]) - sum(I3_entropies[4:6]) + I3_entropies[7]
    return (; subsystem_sizes, S_subsystem_final, I3_partition, I3_entropies, I3_final)
end

"""MPS final EE and I3, without expanding the full state vector.

Two site permutations expose B, C, BC and the disconnected AC as prefixes.
Swaps use `entropy_cutoff` with no maxdim cap, independently of evolution's
bond limit. Normalize the copies after swaps before evaluating bond entropy.
"""
function final_entanglement(model::AnyonModel, state::MPS;
    entropy_cutoff::Float64 = 1e-14,
)
    isfinite(entropy_cutoff) && entropy_cutoff >= 0 ||
        throw(ArgumentError("entropy_cutoff must be finite and nonnegative"))
    L = model.N
    length(state) == L || throw(DimensionMismatch("MPS length must equal model.N"))
    I3_partition = final_entanglement_partition(L)
    A, B, C, D = I3_partition
    normalized = normalize(copy(state))
    subsystem_sizes = collect(1:(L - 1))
    S_subsystem_final = anyon_eelis(model, normalized)
    bcad = normalize(movesites(normalized, vcat(B, C, A, D) .=> collect(1:L);
        cutoff = entropy_cutoff))
    S_B, S_BC = ee_mps(bcad, length(B)), ee_mps(bcad, length(B) + length(C))
    cabd = normalize(movesites(normalized, vcat(C, A, B, D) .=> collect(1:L);
        cutoff = entropy_cutoff))
    S_C, S_AC = ee_mps(cabd, length(C)), ee_mps(cabd, length(C) + length(A))
    I3_entropies = [S_subsystem_final[last(A)], S_B, S_C,
        S_subsystem_final[last(B)], S_AC, S_BC, S_subsystem_final[last(C)]]
    I3_final = sum(I3_entropies[1:3]) - sum(I3_entropies[4:6]) + I3_entropies[7]
    return (; subsystem_sizes, S_subsystem_final, I3_partition, I3_entropies,
        I3_final, entropy_swap_cutoff = entropy_cutoff)
end
