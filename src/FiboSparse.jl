```

    Fibonacci_Ham_sparse(::Type{T}, pbc::Bool=true; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}} -> SparseMatrixCSC{Float64, Int}

# params: a type `T` of BitStr{N, Int}, `pbc` is a boolean value, `measure_class` is a symbol, which default to be :Fibo

# Return the Hamiltonian matrix in anyon basis, which is automatically contain periodic boundary condition (pbc) or open boundary condition (obc).
```
function Fibonacci_Ham_sparse(::Type{T}, pbc::Bool=true;measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}}
    basis=Fibonacci_basis(T,pbc, measure_class=measure_class)

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
Fibonacci_Ham_sparse(N::Int64, pbc::Bool=true) = Fibonacci_Ham_sparse(BitStr{N, Int}, pbc)

```
    Fibonacci_Ham_sparse(::Type{T}, k::Int, Y=nothing; measure_class::Symbol=:Fibo) where {N, T <: BitStr{N}} -> SparseMatrixCSC{ComplexF64, Int}   
    
# params: a type `T` of BitStr{N, Int}, `k` is the momentum of the system, `Y` is the topological charge, which default to be nothing, `measure_class` is a symbol, which default to be :Fibo
# Return the Hamiltonian matrix in given symmetric sector Hilbert space.
```
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
Fibonacci_Ham_sparse(N::Int64, k::Int, Y=nothing) = Fibonacci_Ham_sparse(BitStr{N, Int}, k, Y)

