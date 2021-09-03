using Documenter
using ExactOptimalTransport

makedocs(;
    modules=[ExactOptimalTransport],
    repo="https://github.com/JuliaOptimalTransport/ExactOptimalTransport.jl/blob/{commit}{path}#L{line}",
    sitename="ExactOptimalTransport.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaoptimaltransport.github.io/ExactOptimalTransport.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
    strict=true,
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/JuliaOptimalTransport/ExactOptimalTransport.jl",
    push_preview=true,
    devbranch="main",
)
