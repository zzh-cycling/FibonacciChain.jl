using LinearAlgebra
using FibonacciChain
using JLD
using Arpack
include("./FitEntEntScal.jl")

function Born_trajectories(L::Int64, D::Int64)
    trajlis= Vector{Vector{Int64}}(undef, D)
    for i in 1:D
        trajlis[i] = Int.(rand(Bool, div(L, 2)))
    end
    return trajlis
end

function trajectories_measure(N::Int64, τ::Float64, state::Vector{ET}, trajectories::Vector{Vector{Int64}}) where {ET}
    @assert length(state) == length(Fibonacci_basis(N)) "State vector must have length $(length(Fibonacci_basis(N))), but got $(length(state))"
    for i in eachindex(trajectories)
        traj = trajectories[i]
        if i %2 == 1
            for j in eachindex(traj)
                sign = traj[j] == 0 ? :p : :m
                state=measuremap(N, τ, state, 2j, sign)
                state/=norm(state) # Normalize the state after each measurement
            end
        else
            for j in eachindex(traj)
                sign = traj[j] == 0 ? :p : :m
                state=measuremap(N, τ, state, 2j-1, sign)
                state/=norm(state) # Normalize the state after each measurement
            end
        end
    end
    return state
end