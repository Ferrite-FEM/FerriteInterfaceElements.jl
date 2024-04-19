using TimerOutputs

dto = TimerOutput()
reset_timer!(dto)

using Ferrite
using FerriteInterfaceElements
using Documenter

DocMeta.setdocmeta!(FerriteInterfaceElements, :DocTestSetup, :(using Ferrite, FerriteInterfaceElements); recursive=true)

# Generate tutorials
include("generate.jl")

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
        "Tutorials" => [
            "Tutorials overview" => "tutorials/index.md",
            "tutorials/heat_equation.md",
        ],
        "Reference" => [
            "Reference overview" => "reference/index.md",
            "reference/cells.md",
            "reference/interpolations.md",
            "reference/cellvalues.md",
        ],
        "devdocs/index.md",
    ],
)

deploydocs(;
    repo="github.com/Ferrite-FEM/FerriteInterfaceElements.jl",
    devbranch="main",
    push_preview=true,
)
