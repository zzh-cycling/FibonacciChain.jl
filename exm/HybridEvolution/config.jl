"""
    get_hybrid_time_in_L(p) -> Int

Return the converged exact-dynamics evolution time in units of `L` for the
exclusive hybrid protocol at `γ = 1`. These values were determined from
`L = 12`, 1000-trajectory time-doubling tests of the half-chain entropy,
topological-sector uncertainty, and the full `⟨Y⟩` distribution.
"""
function get_hybrid_time_in_L(p::Real)
    time_in_L = Dict(
        0.1 => 256,
        0.2 => 128,
        0.3 => 64,
        0.4 => 32,
        0.5 => 32,
        0.6 => 32,
        0.7 => 16,
        0.8 => 16,
        0.95 => 8,
    )
    p_float = Float64(p)
    haskey(time_in_L, p_float) || throw(
        ArgumentError(
            "No validated hybrid evolution time for p=$p; " *
            "available probabilities are $(sort!(collect(keys(time_in_L))))",
        ),
    )
    return time_in_L[p_float]
end

"""Return the number of complete hybrid periods for probability `p` and size `L`."""
function get_hybrid_periods(p::Real, L::Integer)
    L >= 1 || throw(ArgumentError("L must be positive, got $L"))
    return get_hybrid_time_in_L(p) * Int(L)
end
