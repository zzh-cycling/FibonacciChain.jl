function Fibonacci_Ham_sparse(::Type{T}, pbc::Bool=true;measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for PXP model, automotically contain pbc or obc
    basis=Fibonacci_basis(T,pbc, measure_class=measure_class)

    l=length(basis)
    # H=spzeros(Float64,(l,l))
    I, J, V = Int[], Int[], Float64[]
    for i in 1:l
        output=actingHam(T, basis[i], pbc)
        states, weights = keys(output), values(output)
        for (idx, m) in enumerate(states)
            j=searchsortedfirst(basis, m)
            push!(I, i); push!(J, j); push!(V, output[m])
        end
    end

    H = sparse(I, J, V, l, l)
    
    return H
end
Fibonacci_Ham_sparse(N::Int64, pbc::Bool=true) = Fibonacci_Ham_sparse(BitStr{N, Int}, pbc)

function Fibonacci_Ham_sparse(::Type{T}, k::Int, Y=nothing; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and topological charge which default to be nothing
#return: the Hamiltonian matrix in given symmetric sector
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 2 or :trivial or :nontrivial, but got $Y"

    basisK, basis_dic = Fibonacci_basis(T, k, Y=Y, measure_class=measure_class)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    # H = spzeros(ComplexF64, (l, l))
    I, J, V = Int[], Int[], ComplexF64[]

    for i in 1:l
        n=basisK[i]
        output = actingHam(T, n, true)
        for m in output
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                # H[i, j] += Yn/Ym * omegak^d
                push!(I, i); push!(J, j); push!(V, Yn/Ym * omegak^d)
            end
        end
    end

    H = sparse(I, J, V, l, l)

    H=(H+H')/2
    if k==0 || k==div(N,2)
        H=real(H)
    end
    return H
end
Fibonacci_Ham_sparse(N::Int64, k::Int, Y=nothing) = Fibonacci_Ham_sparse(BitStr{N, Int}, k, Y)

function iso_total2sec_sparse(::Type{T}, k::Int64, Y=nothing) where {N, T <: BitStr{N}}
#Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 2 or :trivial or :nontrivial, but got $Y"
    @assert k in 0:N-1 "k can only be between from 0 to $(N-1)"
    basis = PXP_basis(T)
    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    # Categorize basis states by their representative
    for i in eachindex(basis)
        state = basis[i]
        category = get_representative(state)[1]
        if haskey(k_dic, category)
            push!(k_dic[category], i)
        else
            k_dic[category] = [i]
        end
    end

    for j in eachindex(basis)
        n=basis[j]
        RS = get_representative(n)[1]
        if RS == n && (k * length(k_dic[RS])) % N == 0
            push!(basisK, n)
        end
    end

    # Initialize sparse matrix
    num_states = length(basis)
    num_categories = length(keys(basisK))
    rows = Vector{Int64}[]
    cols = Vector{Int64}[]
    vals = Vector{Float64}[]

    # Fill the sparse matrix with isometry values
    for (i, state) in enumerate(basisK)
        state_indices = k_dic[state]  
        l = length(state_indices)
        push!(rows, state_indices)
        push!(cols, fill(i, l))
        push!(vals, fill(1.0 / sqrt(l), l))
    end
    
    rows = vcat(rows...)
    cols = vcat(cols...)
    vals = vcat(vals...)
    # Create sparse matrix
    iso_sparse = sparse(rows, cols, vals, num_states, num_categories)

    return iso_sparse
end
iso_total2sec_sparse(N::Int, k::Int64, Y=nothing) = iso_total2sec_sparse(BitStr{N, Int}, k, Y)



# function iso_K2MSS_sparse(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
# #Function to map the MSS basis to the K space basis
#     @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
#     @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

#     basis = PXP_basis(T)
#     basisK, k_dic = PXP_K_basis(T, k)

#     MSS_dic = Dict{Int, Vector{Int64}}()
#     qlist = Vector{Int}(undef, 0)
#     # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
#     if inv==1 && k==0 || k==div(N,2) && inv==-1
#         for i in eachindex(basisK)
#             n = basisK[i]
#             # here we calculate the representative state of the inversion of n
#             nR = get_representative(breflect(n))[1]
#             if n <= min(nR, n)
#                 push!(qlist, length(Set([n, nR])))
#             end
#             n = min(nR, n)
#             if haskey(MSS_dic, n)
#                 push!(MSS_dic[n], i)
#             else
#                 MSS_dic[n] = [i]
#             end
#         end

#     elseif inv==1 && k==div(N,2) || inv==-1 && k==0
#         for i in eachindex(basisK)
#             n = basisK[i]
#             nR = get_representative(breflect(n))[1]
#             if n != nR
#                 n = min(nR, n)
#                 if haskey(MSS_dic, n)
#                     push!(MSS_dic[n], i)
#                 else
#                     MSS_dic[n] = [i]
#                 end
#                 push!(qlist, 2)
#             end
#         end    
        
    
#     end

#     num_states = length(basisK)
#     num_categories = length(keys(MSS_dic))
#     rows = Vector{Int64}[]
#     cols = Vector{Int64}[]
#     vals = Vector{Float64}[]

#     MSS_dic=sort(MSS_dic)
#     for (i, state_indices) in enumerate(values(MSS_dic))
#         l = qlist[i]
#         push!(rows, state_indices)
#         push!(cols, fill(i, l))
#         push!(vals, fill(1.0 / sqrt(l), l))
#     end

#     rows = vcat(rows...)
#     cols = vcat(cols...)
#     vals = vcat(vals...)

#     iso_sparse = sparse(rows, cols, vals, num_states, num_categories)

#     return iso_sparse
# end
# iso_K2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS_sparse(BitStr{N, Int}, k, inv)

# function iso_total2MSS_sparse(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#     # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
#     iso = iso_total2K_sparse(T, k) * iso_K2MSS_sparse(T, k, inv)

#     return iso
# end
# iso_total2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS_sparse(BitStr{N, Int}, k, inv)



