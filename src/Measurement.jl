function measure_basismap(::Type{T}, τ::Float64, state::T, i::Int, sign::Int64, pbc::Bool=true) where {N, T <: BitStr{N}}
    # default for PBC system
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (0, 1) "sign must be either 0 the plus, 1 the minus"
    ϕ = (1+√5)/2
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    
    if τ >= 1e2
        cstτ = 0.5
        coef = sign == 0 ? 0.5 : -0.5
    else
        cstτ = (exp(τ) + 1) / (2 * √(exp(2τ) + 1))
        coef = sign == 0 ? (exp(τ) - 1) / (2 * √(exp(2τ) + 1)) : (1 - exp(τ)) / (2 * √(exp(2τ) + 1))
    end

    if 2<= i <= N-1
        mask=bmask(T,1,2,3) << (N-i-1)
        str100, str101, str010, str001, str000 = T(4) << (N-i-1), T(5) << (N-i-1), T(2) << (N-i-1), T(1) << (N-i-1), T(0) << (N-i-1)
        if state & mask == str000
            return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
        elseif state & mask == str010
            return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
        elseif state & mask == str001
            return state, cstτ+coef
        elseif state & mask == str100
            return state, cstτ+coef
        elseif state & mask == str101
            return state, cstτ-coef
        end
    end
    if pbc
        if i == 1 #count from the left
        mask=bmask(T, N, N-1,1)
        str100, str101, str010, str001, str000 = bmask(T,1), bmask(T, N-1, 1), bmask(T, N), bmask(T, N-1), T(0)
            if state & mask == str000
                return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
            elseif state & mask == str010
                return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
            elseif state & mask == str001
                return state, cstτ+coef
            elseif state & mask == str100
                return state, cstτ+coef
            elseif state & mask == str101
                return state, cstτ-coef
            end
        elseif i == N #count from the left
        mask=bmask(T, N, 2, 1)
        str100, str101, str010, str001, str000 = bmask(T,2), bmask(T, N, 2), bmask(T, 1), bmask(T, N), T(0)
            if state & mask == str000
                return state, X(state,i), cstτ+coef*(1-2ϕ^(-1)), -2*coef*ϕ^(-3/2)
            elseif state & mask == str010
                return state, X(state,i), cstτ+coef*(2ϕ^(-1)-1), -2*coef*ϕ^(-3/2)
            elseif state & mask == str001
                return state, cstτ+coef
            elseif state & mask == str100
                return state, cstτ+coef
            elseif state & mask == str101
                return state, cstτ-coef
            end
        end
    end
end


function measure_matrix(::Type{T}, τ::Float64, idx::Int, sign::Int64, pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    Bmatrix=zeros((l,l))
    for i in 1:l
        outcome = measure_basismap(T, τ, basis[i], idx, sign, pbc)
        if length(outcome) == 4
            outputstate1, outputstate2, output1, output2=outcome
            j2=searchsortedfirst(basis, outputstate2)
            Bmatrix[i,i]+=output1
            Bmatrix[i,j2]+=output2
        else
            outputstate, output=outcome
            Bmatrix[i,i]+=output
        end
    end
    
    return Bmatrix
end

function measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    @assert l == length(state) "state length is expected to be $(l), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        output = measure_basismap(T, τ, basis[i], idx, sign, pbc)
        if length(output) == 4
            outputstate1, outputstate2, output1, output2=output
            j2=searchsortedfirst(basis, outputstate2)
            mapped_state[i]+=output1*state[i] # outputstate1 is the same as basis[i]
            mapped_state[j2]+=output2*state[i]
        else
            outputstate, output1=output # outputstate is the same as basis[i]
            mapped_state[i]+=output1*state[i]
        end
    end
    
    return mapped_state
end
measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true) where {ET} = measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc)

function laddermeasuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    # input a superposition state, and output the braided state
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    mapped_state = zeros(ET, length(state))
    for i in 1:l
        for j in 1:l
            output1 = measure_basismap(T, τ, basis[i], idx, sign, pbc)
            output2 = measure_basismap(T, τ, basis[j], idx, sign, pbc)
            if length(output1) == 4 && length(output2) == 4
                basisi1, basisi2, coefi1, coefi2=output1
                basisj1, basisj2, coefj1, coefj2=output2
                i2=searchsortedfirst(basis, basisi2)
                j2=searchsortedfirst(basis, basisj2)
                # Here noting that the state is a vertorizing density matrix, so the index is i+(j-1)*len, not state[i], state[j]
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi1*coefj2
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj1
                mapped_state[(i2-1)*l+j2]+=state[(i-1)*l+j]*coefi2*coefj2
            elseif length(output1) == 4 && length(output2) == 2
                basisi1, basisi2, coefi1, coefi2=output1
                basisj, coefj=output2
                i2=searchsortedfirst(basis, basisi2)  
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi1*coefj
                mapped_state[(i2-1)*l+j]+=state[(i-1)*l+j]*coefi2*coefj
            elseif length(output1) == 2 && length(output2) == 4
                basisi, coefi=output1
                basisj1, basisj2, coefj1, coefj2=output2
                j2=searchsortedfirst(basis, basisj2)
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj1
                mapped_state[(i-1)*l+j2]+=state[(i-1)*l+j]*coefi*coefj2
            else
                basisi, coefi=output1
                basisj, coefj=output2
                mapped_state[(i-1)*l+j]+=state[(i-1)*l+j]*coefi*coefj
            end
        end
    end
    
    return mapped_state
end
laddermeasuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Int64, pbc::Bool=true) where {ET} = laddermeasuremap(BitStr{N, Int}, τ, state, idx, sign, pbc)

function measurement_enumeration(::Type{T}, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    """
    enumerating all trajectories of measurements on a given initial state.
    
    Args:T, τ, initial_state, measurement_sites, pbc
    
    Returns:
        final_states, trajectories, probabilities
    """
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    
    # Initialize, only one initial state
    current_level_states = [copy(initial_state)]
    current_level_trajectories = [Int64[]]
    current_level_probabilities = [1.0]
    
    for (measurement_idx, site) in enumerate(measurement_sites)
        next_level_states = Vector{Vector{ET}}()
        next_level_trajectories = Vector{Vector{Int64}}()
        next_level_probabilities = Vector{Float64}()
        
        # Branching for each current state
        for (state_idx, state) in enumerate(current_level_states)
            current_trajectory = current_level_trajectories[state_idx]
            current_prob = current_level_probabilities[state_idx]
            
            state_after_p = measuremap(T, τ, state, site, 0, pbc)
            prob_p = state_after_p' * state_after_p  
            
            # if prob_p > 1e-12  
                normalized_state_p = state_after_p / sqrt(prob_p)
                new_trajectory_p = [current_trajectory; 0]
                new_prob_p = current_prob * prob_p
                
                push!(next_level_states, normalized_state_p)
                push!(next_level_trajectories, new_trajectory_p)
                push!(next_level_probabilities, new_prob_p)
            # end
            
            state_after_m = measuremap(T, τ, state, site, 1, pbc)
            prob_m = state_after_m' * state_after_m
            
            # if prob_m > 1e-12 
                normalized_state_m = state_after_m / sqrt(prob_m)
                new_trajectory_m = [current_trajectory; 1]
                new_prob_m = current_prob * prob_m
                
                push!(next_level_states, normalized_state_m)
                push!(next_level_trajectories, new_trajectory_m)
                push!(next_level_probabilities, new_prob_m)
            # end
        end
        
        current_level_states = next_level_states
        current_level_trajectories = next_level_trajectories
        current_level_probabilities = next_level_probabilities
        
        # println("After measurement $(measurement_idx) at site $(site): $(length(current_level_states)) branches")
    end
    

    return current_level_states, current_level_trajectories, current_level_probabilities
end

measurement_enumeration(N::Int, τ::Float64, initial_state::Vector{ET}, measurement_sites::Vector{Int}, pbc::Bool=true) where {ET} = 
measurement_enumeration(BitStr{N, Int}, τ, initial_state, measurement_sites, pbc)


function measurement_tree_visualization(trajectories::Vector{Vector{Int64}}, probabilities::Vector{Float64})
    """
    visualize measurement tree
    """
    total_prob = sum(probabilities)
    normalized_probs = probabilities / total_prob
    
    println("Measurement Tree Visualization:")
    println("==============================")
    
    max_length = maximum(length.(trajectories))
    
    for traj_length in 0:max_length
        level_indices = findall(t -> length(t) == traj_length, trajectories)
        if !isempty(level_indices)
            println("Level $(traj_length):")
            for idx in level_indices
                traj = trajectories[idx]
                prob = normalized_probs[idx]
                indent = "  " ^ traj_length
                if traj_length == 0
                    println("$(indent)Initial state (prob: 1.0)")
                else
                    println("$(indent)$(traj) (prob: $(round(prob, digits=6)))")
                end
            end
            println()
        end
    end
end

function Boundary_measure(::Type{T}, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    
    num_sites = length(measurement_sites)
    
    sample_measured_states = Vector{Vector{Float64}}(undef, num_samples)
    samples = Vector{Vector{Int64}}(undef, num_samples)
    sample_free_energy = Vector{Float64}(undef, num_samples)

    for sample_idx in 1:num_samples
        current_sequence = Vector{Int64}(undef, num_sites)
        current_state = copy(state)  
        total_free_energy = 0.0
        
        # meaure from the left to the right
        for (site_idx, measurement_site) in enumerate(measurement_sites)
           
            state_after_p = measuremap(T, τ, current_state, measurement_site, 0, pbc)
            
            
            prob_p = state_after_p' * state_after_p
            prob_m = 1 - prob_p
            
            random_number = rand(rng)
            if random_number < prob_p
                current_sequence[site_idx] = 0
                current_state = state_after_p ./ sqrt(prob_p)
                total_free_energy += -log(prob_p)
            else
                state_after_m = measuremap(T, τ, current_state, measurement_site, 1, pbc)
                current_sequence[site_idx] = 1
                current_state = state_after_m ./ sqrt(prob_m)
                total_free_energy += -log(prob_m)
            end
        end
        
        sample_measured_states[sample_idx] = current_state
        samples[sample_idx] = current_sequence
        sample_free_energy[sample_idx] = total_free_energy
    end
    
    return sample_measured_states, samples, sample_free_energy
end

Boundary_measure(N::Int, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int},num_samples::Int=1000, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true) where {ET} = Boundary_measure(BitStr{N, Int}, τ, state, measurement_sites, num_samples, rng, pbc)

function Boundarypost_selection(N::Int64, τ::Float64, state::Vector{ET}, measurement_sites::Vector{Int}, sign::Int64, pbc::Bool=true) where {ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    num_sites = length(measurement_sites)

    current_sequence = Vector{Int64}(undef, num_sites)
    current_state = copy(state)  
    total_free_energy = 0.0
    
    # meaure from the left to the right
    for (site_idx, measurement_site) in enumerate(measurement_sites)
       
        state_after_measure = measuremap(N, τ, current_state, measurement_site, sign, pbc)

        prob = state_after_measure' * state_after_measure

        current_sequence[site_idx] = sign
        current_state = state_after_measure ./ sqrt(prob)
        total_free_energy += -log(prob)
    end
    
    return current_state, current_sequence, total_free_energy
end

function Bulkmeasure(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, rng::MersenneTwister=MersenneTwister(), pbc::Bool=true) where {ET}
    @assert length(state) == length(Fibonacci_basis(N)) "State vector must have length $(length(Fibonacci_basis(N))), but got $(length(state))"
    sample = zeros(Int, D, div(N,2))
    sample_free_energy = Vector{Float64}(undef, D)
    sample_measured_states = Vector{Vector{ET}}(undef, D)

    current_state = copy(state)  
    
    for layer in 1:D
        current_sequence = zeros(Int, div(N,2))
        total_free_energy = 0.0
        
        # Alternating measurement pattern
        if layer % 2 == 1
            measurement_sites = collect(2:2:N)  # odd layers: even sites
        else
            measurement_sites = collect(1:2:N)  # even layers: odd sites
        end
        
        measurement_τ = (layer == D) ? τ/2 : τ
        # measure :sqrtp is Pi 0/tau, measure :sqrtm is Pi 1, measure 0 is Pi 0/tau, measure 1 is Pi 1
            
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            state_after_p = measuremap(N, measurement_τ, current_state, measurement_site, 0, pbc)
            
            prob_sqrtp = state_after_p' * state_after_p
            prob_sqrtm = 1 - prob_sqrtp
            random_number = rand(rng)  
            if random_number < prob_sqrtp
                current_sequence[site_idx] = 0
                current_state = state_after_p ./ sqrt(prob_sqrtp)
                total_free_energy += -log(prob_sqrtp)
            else
                state_after_m = measuremap(N, measurement_τ, current_state, measurement_site, 1, pbc)
                current_sequence[site_idx] = 1
                current_state = state_after_m ./ sqrt(prob_sqrtm)
                total_free_energy += -log(prob_sqrtm)
            end
        end

        sample_measured_states[layer] = current_state
        sample[layer, :] = current_sequence
        sample_free_energy[layer] = total_free_energy
    end

    return sample_measured_states, sample, sample_free_energy
end

function Bulkpost_selection(N::Int64, τ::Float64, state::Vector{ET}, D::Int64, sign::Int64, pbc::Bool=true) where {ET}
    @assert length(state) == length(Fibonacci_basis(N)) "State vector must have length $(length(Fibonacci_basis(N))), but got $(length(state))"
    # N is the number of sites, τ is the measurement parameter, state is the initial state vector, D is the layer depth of the measurement tree
    sample = Vector{Vector{Int64}}(undef, D)
    sample_free_energy = Vector{Float64}(undef, D)
    sample_measured_states = Vector{Vector{ET}}(undef, D)

    current_state = copy(state)  

    for layer in 1:D
        current_sequence = Vector{Int64}(undef, div(N,2))
        total_free_energy = 0.0
        # total_free_energy is the log probability of the sample (average free energy), so it should be initialized to 0.0
        
        # Alternating measurement pattern
        if layer % 2 == 1
            measurement_sites = collect(2:2:N)  # odd layers: even sites
        else
            measurement_sites = collect(1:2:N)  # even layers: odd sites
        end
        
        measurement_τ = (layer == D) ? τ/2 : τ
        
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            state_after_p = measuremap(N, measurement_τ, current_state, measurement_site, sign, pbc)
            current_sequence[site_idx] = sign
            prob_p = state_after_p' * state_after_p
            current_state = state_after_p ./ sqrt(prob_p)
            total_free_energy += -log(prob_p)
        end

        sample_measured_states[layer] = current_state
        sample[layer] = current_sequence
        sample_free_energy[layer] = total_free_energy
        
    end

    return sample_measured_states, sample, sample_free_energy
end

# Helper function to apply measurements to a layer
function apply_measurement_layer!(N::Int64, state::Vector{T}, τ::Float64, layer_sample::Vector{Int64}, layer_idx::Int64, pbc::Bool=true) where {T}
    if layer_idx % 2 == 1
        measurement_sites = collect(2:2:N)  # odd sites anyons, even sites qubits
    else
        measurement_sites = collect(1:2:N)  # even sites anyons, odd sites qubits
    end
    for (idx, measurement_type) in enumerate(layer_sample)
        state = measuremap(N, τ, state, measurement_sites[idx], measurement_type, pbc)
        normalize!(state)
    end
    return state
end

function generate_state(τ::Float64, state::Vector{T}, sample::ET, temp::Bool=false, pbc::Bool=true) where{T, ET}

    if ET == Vector{Int}
        N = 2 * length(sample)
        return apply_measurement_layer!(N, state, τ, sample, 1, pbc)
        
    elseif ET == Matrix{Int}
        D, N = size(sample, 1), 2 * size(sample, 2)
        statelis = temp ? Vector{Vector{T}}(undef, D) : nothing
        # if ET is Vector{Int64} and temp is true, we return temporary states.
        for layer in 1:D
            τ_eff = (layer == D) ? τ/2 : τ
            state = apply_measurement_layer!(N, state, τ_eff, sample[layer, :], layer, pbc)
            
            if temp
                statelis[layer] = copy(state)
            end
        end
        
        return temp ? statelis : state
    end
end