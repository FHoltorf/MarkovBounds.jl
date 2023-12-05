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
    warnonly = Documenter.except(),
    pages=[
        "Home" => "index.md",
        "Background" => "background.md",
        "Tutorials" => ["tutorials/birth_death_process.md",
                        "tutorials/jump_diffusion_process.md",
                        "tutorials/optimal_control.md"],
        "Functionalities" => "functionalities.md" 
    ],
)

deploydocs(;
    repo="github.com/FHoltorf/MarkovBounds.jl",
)
