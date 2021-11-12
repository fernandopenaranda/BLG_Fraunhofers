using BLG_Fraunhofers
using Documenter

DocMeta.setdocmeta!(BLG_Fraunhofers, :DocTestSetup, :(using BLG_Fraunhofers); recursive=true)

makedocs(;
    modules=[BLG_Fraunhofers],
    authors="Fernando Peñaranda",
    repo="https://github.com/fernandopenaranda/BLG_Fraunhofers.jl/blob/{commit}{path}#{line}",
    sitename="BLG_Fraunhofers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fernandopenaranda.github.io/BLG_Fraunhofers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/BLG_Fraunhofers.jl",
    devbranch="main",
)
