using FibonacciChain
using JLD2
using LinearAlgebra

"""
    kw_postselection(L, periods; tau=log(1+sqrt(2)), outcome=false, half_boundary=true)

Evolve the all-plus state under a fixed-outcome Ising measurement record.
Returns KW expectations and the final eigenstate residual. `half_boundary`
uses the package's final sqrt-X temporal boundary; earlier points are bulk periods.
No files are written by this function.
"""
function kw_postselection(L::Integer, periods::Integer;
    tau::Real=log(1+sqrt(2)), outcome::Bool=false, half_boundary::Bool=true)
    L >= 2 || throw(ArgumentError("L must be at least 2"))
    periods >= 1 || throw(ArgumentError("periods must be positive"))
    isfinite(tau) && tau > 0 || throw(ArgumentError("tau must be positive and finite"))
    model = AnyonModel(SpinHalf(), Int(L); model_type=:Ising, pbc=true, measure_operator=:X)
    state = fill(2.0^(-L/2), 2^L)
    record = BitMatrix(fill(outcome, 2periods, FibonacciChain._samples_per_layer(model)))
    config = MeasureConfig(τ=Float64(tau), mode=:sample, t₂=Int(periods),
        track_y_expectation=true, enable_τ_eff=half_boundary)
    result = bulk_evolution(model, state, config, record)
    v = result.state
    Dv = kramers_wannier_map(model) * v
    expectation = dot(v, Dv) / dot(v, v)
    residual = norm(Dv - expectation*v) / norm(v)
    return (; L=Int(L), periods=Int(periods), τ=Float64(tau), gamma=tanh(tau),
        initial_state="all_plus", record=outcome ? "all_ones" : "all_zeros",
        kw_expectation_values=Float32.(result.y_expectation_values),
        final_kw_expectation=expectation, final_eigenstate_residual=residual,
        boundary=half_boundary ? "sqrt-X half layer on the last layer (enable_τ_eff=true)" : "full layers")
end

"""Write one explicit output file; refuse to overwrite existing scientific data."""
function save_kw_postselection(output::AbstractString, L::Integer, periods::Integer; kwargs...)
    ispath(output) && throw(ArgumentError("Output already exists: $output"))
    result = kw_postselection(L, periods; kwargs...)
    mkpath(dirname(abspath(output)))
    jldsave(output; result...)
    return output
end

# From the FibonacciChain root:
# julia --project=. exm/Born_Ising/kw_postselection.jl OUTPUT.jld2 L PERIODS [TAU]
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (3,4) || error("Usage: kw_postselection.jl OUTPUT.jld2 L PERIODS [TAU]")
    tau = length(ARGS) == 4 ? parse(Float64, ARGS[4]) : log(1+sqrt(2))
    println(save_kw_postselection(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]); tau))
end
