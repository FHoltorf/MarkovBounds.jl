using MarkovBoundsSOS
using Documenter

DocMeta.setdocmeta!(MarkovBoundsSOS, :DocTestSetup, :(using MarkovBoundsSOS); recursive=true)

makedocs(;
    modules=[MarkovBoundsSOS],
    authors="Flemming Holtorf",
    repo="https://github.com/FHoltorf/MarkovBoundsSOS.jl/blob/{commit}{path}#{line}",
    sitename="MarkovBoundsSOS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://FHoltorf.github.io/MarkovBoundsSOS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FHoltorf/MarkovBoundsSOS.jl",
)
