module FibonacciChain

using BitBasis, LinearAlgebra

export Fibonacci_Ham, Fibonacci_ferroHam, Fibonacci_basis, rdm_Fibo, ladderrdm
export eelis_Fibo_state

include("basis.jl")
include("Observable.jl")
end
