"""
    anyon_ham_sparse(model::AnyonModel; kwargs...)

Construct the anyon chain Hamiltonian as a sparse matrix.

# Arguments
- `model::AnyonModel`: An `AnyonModel` object containing the system parameters.
- `kwargs...`: Additional keyword arguments passed to `actingHam`, e.g., `J`, `h` for the Ising model.

# Returns
- `SparseMatrixCSC{Float64, Int}`: The sparse Hamiltonian matrix in the anyon basis.

# Examples
```jldoctest
julia> using FibonacciChain, BitBasis, SparseArrays, LinearAlgebra

julia> N = 6;
julia> model = AnyonModel(FibonacciAnyon(), N, pbc=true);
julia> H_sparse = anyon_ham_sparse(model);

julia> H_sparse isa SparseMatrixCSC # Check it's a sparse matrix
true

julia> ishermitian(H_sparse) # Should be Hermitian
true

julia> H_dense = anyon_ham(model); # Compare with dense version for small system

julia> norm(Matrix(H_sparse) - H_dense) < 1e-10
true
```
"""
function anyon_ham_sparse(model::AnyonModel; kwargs...)
    basis=anyon_basis(model)

    l=length(basis)
    I, J, V = Int[], Int[], Float64[]
    for i in 1:l
        output=actingHam(model, basis[i]; kwargs...)
        states, weights = keys(output), values(output)
        for (idx, m) in enumerate(states)
            j=searchsortedfirst(basis, m)
            push!(I, i); push!(J, j); push!(V, output[m])
        end
    end

    H = sparse(I, J, V, l, l)
    
    return H
end


"""
    anyon_ham_sparse(model::AnyonModel, k::Int; symmetric_block=nothing, kwargs...)

Construct the Hamiltonian in a specific momentum and/or topological sector as a sparse matrix.

# Arguments
- `model::AnyonModel`: An `AnyonModel` object containing the system parameters.
- `k::Int`: The momentum quantum number (0 ≤ k ≤ N-1).
- `symmetric_block`: Optional. The topological charge sector (e.g., `:trivial`, `:nontrivial`).
- `kwargs...`: Additional keyword arguments passed to `actingHam`, e.g., `J`, `h`.

# Returns
- `SparseMatrixCSC`: The sparse Hamiltonian in the specified symmetric sector. The element type will be `ComplexF64` for general `k` and `Float64` for `k=0` or `k=N/2`.
"""
function anyon_ham_sparse(model::AnyonModel, k::Int; symmetry_block=nothing, kwargs...)
#params: a int of lattice number, momentum of system and topological charge which default to be nothing
#return: the Hamiltonian matrix in given symmetric sector
    N = model.N
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"
    @assert symmetry_block === nothing || symmetry_block in [0, 1, :tau, :trivial] "Y is expected to be nothing or 1 or 2 or :trivial or :nontrivial, but got $Y"

    basisK, basis_dic = anyon_basis(model, k; symmetry_block= symmetry_block)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    # H = spzeros(ComplexF64, (l, l))
    I, J, V = Int[], Int[], ComplexF64[]

    for i in 1:l
        n=basisK[i]
        output = actingHam(model, n; kwargs...)
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

    if k==0 || k==div(N,2)
        H=real(H)
    end
    H=(H+H')/2
    return H
end
