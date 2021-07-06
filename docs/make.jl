using MarkovBounds
using Documenter

DocMeta.setdocmeta!(MarkovBounds, :DocTestSetup, :(using MarkovBounds); recursive=true)

makedocs(;
    modules=[MarkovBounds],
    authors="Flemming Holtorf",
    repo="https://github.com/FHoltorf/MarkovBounds.jl/blob/{commit}{path}#{line}",
    sitename="MarkovBounds.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://FHoltorf.github.io/MarkovBounds.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FHoltorf/MarkovBounds.jl",
)
