using Ferrite, FerriteInterfaceElements, OrderedCollections

@testset "VTK Export" begin
    mktempdir() do tmp
        for grid in [
            generate_grid(Triangle, (4,4)),
            generate_grid(Quadrilateral, (4,4)),
            generate_grid(Tetrahedron, (4,4,4)),
            generate_grid(Hexahedron, (4,4,4)),
        ]
            addcellset!(grid, "A", x->norm(x,Inf) ≤ 0.5)
            addcellset!(grid, "B", setdiff(OrderedSet(1:getncells(grid)), getcellset(grid, "A")))
            grid2 = insert_interfaces(grid, ["A", "B"])
            VTKGridFile("debug.vtu", grid2) do vtk
                Ferrite.write_cellset(vtk, grid2)
            end
        end
    end
end
