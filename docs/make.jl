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
        "Manual" => [
            "Basis Functions" => "basis.md",
            "Observables" => "observables.md", 
            "Measurements" => "measurements.md",
            "MPS Methods" => "mps.md",
            "Examples" => "examples.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/zzh-cycling/FibonacciChain.jl",
    devbranch="main",
)
