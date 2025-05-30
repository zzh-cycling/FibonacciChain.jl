function ee(subrm::Matrix{ET}) where {ET}
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
    @assert ishermitian(subrm) "The reduced density matrix is not hermitian."
    spectrum=eigvals(subrm)
    EE=0
    for i in eachindex(spectrum)
        v=abs(spectrum[i])
            if v>1e-8
                EE+=-v*log(v)
            end
    end

    return EE
end

function eelis_Fibo_state(N::Int64,splitlis::Vector{Int64},state::Vector{ET},pbc::Bool=true) where {ET}
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_Fibo(N, collect(1:splitlis[m]), state, pbc)
        EE_lis[m]=ee(subrho)
    end
    return EE_lis
end


