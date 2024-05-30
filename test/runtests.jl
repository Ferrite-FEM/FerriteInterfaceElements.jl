using FerriteInterfaceElements
using FerriteInterfaceElements.Ferrite
using Test

@testset "FerriteInterfaceElements.jl" begin
    include("test_cells.jl")
    
    include("test_interpolations.jl")

    include("test_cellvalues.jl")

    include("test_grid.jl")

    include("test_examples.jl")
end
