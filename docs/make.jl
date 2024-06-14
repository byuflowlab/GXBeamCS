using GXBeamCS
using Documenter

DocMeta.setdocmeta!(GXBeamCS, :DocTestSetup, :(using GXBeamCS); recursive=true)

makedocs(;
    modules=[GXBeamCS],
    authors="Andrew Ning <aning@byu.edu>, Adam Cardoza",
    sitename="GXBeamCS.jl",
    format=Documenter.HTML(;
        canonical="https://byuflowlab.github.io/GXBeamCS.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/byuflowlab/GXBeamCS.jl",
    devbranch="main",
)
