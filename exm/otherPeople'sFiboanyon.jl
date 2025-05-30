using LinearAlgebra

function Chains(n::Int)
    n <= 0 && return Matrix{Int}(undef, 0, 0)
    chains = [[0], [1]]
    n == 1 && return hcat(chains...)'

    output = Vector{Vector{Int}}()
    for _ in 1:(n-1)
        empty!(output)
        for chain in chains
            if chain[end] == 1
                push!(output, vcat(chain, 0))
                push!(output, vcat(chain, 1))
            else
                push!(output, vcat(chain, 1))
            end
        end
        chains = copy(output)
    end
    return permutedims(hcat(chains...))
end



function Periodic_Chains(n::Int)
    pc = Vector{Vector{Int}}()
    c = Chains(n + 1)
    for i in eachrow(c)
        if i[1] == i[end]
            push!(pc, i[1:end-1])
        end
    end
    return permutedims(hcat(pc...))
end

function Null(m::Int)
    return zeros(Int, m, m)
end

function RRotation(Space1::Matrix{Int})
    space2 = Vector{Vector{Int}}()
    for a in eachrow(Space1)
        b = circshift(a, -1)
        push!(space2, b)
    end
    return permutedims(hcat(space2...))
end

function LRotation(Space1::Matrix{Int})
    space2 = Vector{Vector{Int}}()
    for a in eachrow(Space1)
        b = circshift(a, 1)
        push!(space2, b)
    end
    return permutedims(hcat(space2...))
end

function Comp_Chains(n::Int, k::Int)
    a = Vector{Vector{Int}}()
    c = Chains(n - 1)
    for _ in 1:k
        c = RRotation(c)
    end
    for i in eachrow(c)
        if i[1] != i[end]
            push!(a, vcat([1], i))
        else
            push!(a, vcat([0], i))
            if i[1] == 1
                push!(a, vcat([1], i))
            end
        end
    end
    for _ in 1:k
        a = LRotation(permutedims(hcat(a...)))
    end
    return permutedims(hcat(a...))
end

function Reorder(Space1::Matrix{Int}, Space2::Matrix{Int})
    c = Vector{Vector{Int}}()
    spec = 3000
    for (i, a) in enumerate(eachrow(Space1))
        p = 0
        for b in eachrow(Space2)
            if a == b
                push!(c, a)
                p = 1
                break
            end
        end
        if p == 0
            spec = i
        end
    end
    if spec != 3000
        for b in eachrow(Space2)
            p = 0
            for a in eachrow(Space1)
                if a == b
                    p = 1
                    break
                end
            end
            if p == 0
                insert!(c, spec, b)
            end
        end
    end
    return permutedims(hcat(c...))
end

function Fusion(n::Int)
    # Placeholder for F_YM; replace with actual definition
    F_YM = [1.0 0.0; 0.0 1.0]  # Example 2x2 matrix
    space2 = Chains(n - 1)
    space1 = Periodic_Chains(n)
    H = zeros(ComplexF64, size(space1, 1), size(space1, 1))
    for i in 0:n-1
        F = zeros(ComplexF64, size(space1, 1), size(space1, 1))
        P = zeros(ComplexF64, size(space1, 1), size(space1, 1))
        space3 = Vector{Vector{Int}}()
        for a_1 in eachrow(Periodic_Chains(n))
            a_2 = circshift(a_1, -i)
            push!(space3, a_2)
        end
        space3 = permutedims(hcat(space3...))
        space4 = Reorder(Periodic_Chains(n), Comp_Chains(n, i))
        for c in eachrow(space2)
            space1 = Periodic_Chains(n)
            if c[1] != c[end]
                for (k, d) in enumerate(eachrow(space3))
                    if c == d[2:end]
                        space1[k, i+1] = 1
                        for (g_ind, g) in enumerate(eachrow(space4))
                            if space1[k, :] == g
                                F[k, g_ind] = 1
                                P[g_ind, g_ind] = -1
                            end
                        end
                    end
                end
                elseif c[1] == 1 && c[end] == 1
                    for (k, d) in enumerate(eachrow(space3))
                        if c == d[2:end]
                            space1[k, i+1] = 0
                            for (g_ind, g) in enumerate(eachrow(space4))
                                if space1[k, :] == g
                                    F[k, g_ind] = F_YM[d[1]+1, 1]
                                end
                            end
                            space1[k, i+1] = 1
                            for (g_ind, g) in enumerate(eachrow(space4))
                                if space1[k, :] == g
                                    F[k, g_ind] = F_YM[d[1]+1, 2]
                                    P[g_ind, g_ind] = -1
                                end
                            end
                        end
                    end
                
                    elseif c[1] == 0 && c[end] == 0
                        for (k, d) in enumerate(eachrow(space3))
                            if c == d[2:end]
                                space1[k, i+1] = 0
                                for (g_ind, g) in enumerate(eachrow(space4))
                                    if space1[k, :] == g
                                        F[k, g_ind] = 1
                                    end
                                end
                            end
                        end
                    end
                    println(F)
                    H += inv(F) * P * F
                end
                return H
            end
        end