using AutoChem
using Documenter

DocMeta.setdocmeta!(AutoChem, :DocTestSetup, :(using AutoChem); recursive=true)

makedocs(;
    modules=[AutoChem],
    authors="John Waczak <john.louis.waczak@gmail.com>",
    repo="https://github.com/john-waczak/AutoChem.jl/blob/{commit}{path}#{line}",
    sitename="AutoChem.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://john-waczak.github.io/AutoChem.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/john-waczak/AutoChem.jl",
    devbranch="main",
)
