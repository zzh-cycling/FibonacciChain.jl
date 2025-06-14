function measure_basismap(::Type{T}, τ::Float64, state::T, i::Int, sign::Symbol, pbc::Bool=true) where {N, T <: BitStr{N}}
    # default for PBC system
    @assert 1 <= i <= N "Index i must be in the range [1, N]"
    @assert sign in (:p, :m) "sign must be either :p the plus or :m the minus"
    ϕ = (1+√5)/2
    fl=bmask(T, N)
    X(state,i) = flip(state, fl >> (i-1))
    cstτ = (exp(τ)+1)/2√(exp(2τ)+1)
    if sign == :p
        coef = (exp(τ)-1)/2√(exp(2τ)+1)
    else
        coef = (1-exp(τ))/2√(exp(2τ)+1)
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


function measure_matrix(::Type{T}, τ::Float64, idx::Int, sign::Symbol, pbc::Bool=true) where {N, T <: BitStr{N}}
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

function measuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
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
measuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {ET} = measuremap(BitStr{N, Int}, τ, state, idx, sign, pbc)

function laddermeasuremap(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
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
laddermeasuremap(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {ET} = laddermeasuremap(BitStr{N, Int}, τ, state, idx, sign, pbc)

function Sampling(::Type{T}, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    @assert pbc || (2 <= idx <= N-1) "Index idx must be in the range [2, N-1] for open boundary conditions"
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"

    basis=Fibonacci_basis(T, pbc)
    l=length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    mapped_state = zeros(ET, length(state))

    Problis=zeros(length(2:2:N))
    for i in 2:2:N
        Stp = laddermeasuremap(T, τ, state, i, :p, pbc)
        Stm = laddermeasuremap(T, τ, state, i, :m, pbc)
    end
end
Sampling(N::Int, τ::Float64, state::Vector{ET}, idx::Int, sign::Symbol, pbc::Bool=true) where {ET} = Sampling(BitStr{N, Int}, τ, state, idx, sign, pbc)

function Sampling(::Type{T}, τ::Float64, state::Vector{ET}, num_samples::Int=1000, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    @assert ET != Int "The state should be a Float or Complex list, not an integer list"
    
    basis = Fibonacci_basis(T, pbc)
    l = length(basis)
    @assert l^2 == length(state) "state length is expected to be $(l^2), but got $(length(state))"
    
    measurement_sites = collect(2:2:N) 
    num_sites = length(measurement_sites)
    
    samples = Vector{Vector{Symbol}}(undef, num_samples)
    sample_weights = Vector{Float64}(undef, num_samples)
    
    for sample_idx in 1:num_samples
        # 初始化当前采样序列和当前状态
        current_sequence = Vector{Symbol}(undef, num_sites)
        current_state = copy(state)  # 从输入的密度矩阵态开始
        total_weight = 1.0
        
        # 从左到右顺序采样每个测量位置
        for (site_idx, measurement_site) in enumerate(measurement_sites)
            # 计算在当前状态下，该位置测量 :p 和 :m 的概率
            
            # 应用 :p 测量后的状态
            state_after_p = laddermeasuremap(T, τ, current_state, measurement_site, :p, pbc)
            # 应用 :m 测量后的状态  
            state_after_m = laddermeasuremap(T, τ, current_state, measurement_site, :m, pbc)
            
            # 计算概率权重（trace norm）
            prob_p = real(sum(state_after_p))  # Tr(ρ_p)
            prob_m = real(sum(state_after_m))  # Tr(ρ_m)
            
            # 归一化概率
            total_prob = prob_p + prob_m
            if total_prob ≈ 0
                # 如果总概率接近0，随机选择
                prob_p_normalized = 0.5
            else
                prob_p_normalized = prob_p / total_prob
            end
            
            # 随机采样决定测量结果
            random_number = rand()
            if random_number < prob_p_normalized
                # 选择 :p 结果
                current_sequence[site_idx] = :p
                current_state = state_after_p / prob_p  # 归一化状态
                total_weight *= prob_p
            else
                # 选择 :m 结果
                current_sequence[site_idx] = :m
                current_state = state_after_m / prob_m  # 归一化状态
                total_weight *= prob_m
            end
        end
        
        # 存储采样结果
        samples[sample_idx] = current_sequence
        sample_weights[sample_idx] = total_weight
    end
    
    return samples, sample_weights
end

# 便利函数
Sampling(N::Int, τ::Float64, state::Vector{ET}, num_samples::Int=1000, pbc::Bool=true) where {ET} = 
    Sampling(BitStr{N, Int}, τ, state, num_samples, pbc)

# 计算采样后的物理量平均值
function compute_sample_average(samples::Vector{Vector{Symbol}}, 
                               sample_weights::Vector{Float64},
                               observable_func::Function)
    """
    计算物理量在采样序列上的Born轨迹平均值
    
    Args:
        samples: 采样得到的测量序列
        sample_weights: 对应的权重（虽然在重要采样中通常不需要显式使用）
        observable_func: 计算物理量的函数，输入为测量序列，输出为物理量值
    
    Returns:
        平均值
    """
    num_samples = length(samples)
    total = 0.0
    
    for i in 1:num_samples
        observable_value = observable_func(samples[i])
        total += observable_value
    end
    
    return total / num_samples
end

# 示例：计算纠缠熵的采样平均
function entanglement_entropy_from_sequence(::Type{T}, τ::Float64, sequence::Vector{Symbol}, 
                                          initial_state::Vector{ET}, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    """
    根据测量序列计算最终态的纠缠熵
    """
    measurement_sites = collect(2:2:N)
    current_state = copy(initial_state)
    
    # 按序列应用所有测量
    for (site_idx, measurement_site) in enumerate(measurement_sites)
        sign = sequence[site_idx]
        current_state = laddermeasuremap(T, τ, current_state, measurement_site, sign, pbc)
        # 归一化
        norm = real(sum(current_state))
        if norm > 1e-12
            current_state /= norm
        end
    end
    
    # 计算纠缠熵（这里需要根据具体的纠缠熵计算方法实现）
    # 例如，可以重塑为密度矩阵并计算von Neumann熵
    basis = Fibonacci_basis(T, pbc)
    l = length(basis)
    ρ = reshape(current_state, l, l)
    
    # 计算约化密度矩阵的纠缠熵（这里简化处理）
    eigenvals = real.(eigvals(ρ))
    eigenvals = eigenvals[eigenvals .> 1e-12]  # 去除数值噪声
    
    if length(eigenvals) == 0
        return 0.0
    end
    
    # von Neumann 熵
    entropy = -sum(eigenvals .* log.(eigenvals))
    return entropy
end

# 使用示例
function run_importance_sampling_analysis(::Type{T}, τ::Float64, initial_state::Vector{ET}, 
                                        num_samples::Int=10000, pbc::Bool=true) where {N, T <: BitStr{N}, ET}
    """
    运行重要采样分析，计算纠缠熵的Born轨迹平均
    """
    # 执行重要采样
    samples, weights = Sampling(T, τ, initial_state, num_samples, pbc)
    
    # 定义纠缠熵计算函数
    entropy_func = (seq) -> entanglement_entropy_from_sequence(T, τ, seq, initial_state, pbc)
    
    # 计算平均纠缠熵
    avg_entropy = compute_sample_average(samples, weights, entropy_func)
    
    println("Average entanglement entropy from $(num_samples) samples: $(avg_entropy)")
    
    # 计算其他统计量
    sequence_frequencies = Dict{Vector{Symbol}, Int}()
    for seq in samples
        sequence_frequencies[seq] = get(sequence_frequencies, seq, 0) + 1
    end
    
    println("Number of unique sequences sampled: $(length(sequence_frequencies))")
    println("Most frequent sequences:")
    sorted_freqs = sort(collect(sequence_frequencies), by=x->x[2], rev=true)
    for i in 1:min(5, length(sorted_freqs))
        seq, freq = sorted_freqs[i]
        println("  $(seq): $(freq) times ($(freq/num_samples*100)%)")
    end
    
    return samples, weights, avg_entropy
end