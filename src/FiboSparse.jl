"""
    anyon_ham_sparse(::Type{T}, pbc::Bool=true; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Construct Fibonacci chain Hamiltonian as sparse matrix.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `pbc::Bool=true`: Periodic boundary conditions
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `SparseMatrixCSC{Float64, Int}`: Hamiltonian matrix in anyon basis
"""
function anyon_ham_sparse(::Type{T}, pbc::Bool=true;anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T,pbc, anyon_type=anyon_type)

    l=length(basis)
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
anyon_ham_sparse(N::Int64, pbc::Bool=true) = anyon_ham_sparse(BitStr{N, Int}, pbc)

"""
    anyon_ham_sparse(::Type{T}, k::Int, Y=nothing; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}

Construct Hamiltonian in momentum-topological sector as sparse matrix.

# Arguments
- `T::Type`: BitStr type specifying chain length N
- `k::Int`: Momentum quantum number (0 ≤ k ≤ N-1)
- `Y`: Topological charge sector
- `anyon_type::Symbol=:Fibo`: Model type

# Returns
- `SparseMatrixCSC{ComplexF64, Int}`: Hamiltonian in symmetric sector
"""
function anyon_ham_sparse(::Type{T}, k::Int, Y=nothing; anyon_type::Symbol=:Fibo) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and topological charge which default to be nothing
#return: the Hamiltonian matrix in given symmetric sector
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert Y === nothing || Y in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 2 or :trivial or :nontrivial, but got $Y"

    basisK, basis_dic = Fibonacci_basis(T, k, Y=Y, anyon_type=anyon_type)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    # H = spzeros(ComplexF64, (l, l))
    I, J, V = Int[], Int[], ComplexF64[]

    for i in 1:l
        n=basisK[i]
        output = actingHam(T, n, true)
        states, weights = keys(output), values(output)
        for m in states
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                # H[i, j] += Yn/Ym * omegak^d
                push!(I, i); push!(J, j); push!(V, Yn/Ym * omegak^d*output[m])
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
anyon_ham_sparse(N::Int64, k::Int, Y=nothing) = anyon_ham_sparse(BitStr{N, Int}, k, Y)

