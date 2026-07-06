using Ferrite, FerriteInterfaceElements, OrderedCollections

getinterfaceshape(cell::Triangle) = RefLine
getinterfaceshape(cell::Quadrilateral) = RefLine
getinterfaceshape(cell::Tetrahedron) = RefTriangle
getinterfaceshape(cell::Hexahedron) = RefQuadrilateral

@testset "VTK Export Smoke Test" begin
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
            dh = DofHandler(grid2)
            set_bulk = union(getcellset(grid2, "A"), getcellset(grid2, "B"))
            add!(SubDofHandler(dh, set_bulk),      :u, Lagrange{Ferrite.getrefshape(grid.cells[1]), 1}())
            set_interface = getcellset(grid2, "interfaces")
            add!(SubDofHandler(dh, set_interface), :u, InterfaceCellInterpolation(Lagrange{getinterfaceshape(grid.cells[1]), 1}()))
            close!(dh)
            VTKGridFile("debug.vtu", grid2) do vtk
                Ferrite.write_solution(vtk, dh, rand(ndofs(dh)))
                Ferrite.write_cellset(vtk, grid2)
            end
        end
    end
end
