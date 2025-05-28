using FibonacciChain
using Documenter

DocMeta.setdocmeta!(FibonacciChain, :DocTestSetup, :(using FibonacciChain); recursive=true)

makedocs(;
    modules=[FibonacciChain],
    authors="Zhaohui Zhi",
    sitename="FibonacciChain.jl",
    format=Documenter.HTML(;
        canonical="https://zzh-cycling.github.io/FibonacciChain.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zzh-cycling/FibonacciChain.jl",
    devbranch="main",
)
