using FerriteInterfaceElements
using Documenter

DocMeta.setdocmeta!(FerriteInterfaceElements, :DocTestSetup, :(using FerriteInterfaceElements); recursive=true)

makedocs(;
    modules=[FerriteInterfaceElements],
    authors="Kim Louisa Auth <kim.auth@chalmers.se>, Elias Börjesson <elias.borjesson@chalmers.se>, David Rollin <d.rollin@tu-braunschweig.de>",
    sitename="FerriteInterfaceElements.jl",
    format=Documenter.HTML(;
        canonical="https://Ferrite-FEM.github.io/FerriteInterfaceElements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Ferrite-FEM/FerriteInterfaceElements.jl",
    devbranch="main",
)
