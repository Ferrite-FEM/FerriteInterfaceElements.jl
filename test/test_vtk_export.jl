using Ferrite, FerriteInterfaceElements, OrderedCollections

@testset "VTK export" begin
    @testset "InterfaceCell to VTK cell" begin
        # 2D: the VTK quad is in cyclic order, the interface cell is here-then-there
        cell = InterfaceCell(Line((1, 2)), Line((4, 3)))
        @test Ferrite.cell_to_vtkcell(typeof(cell)) == Ferrite.VTKCellTypes.VTK_QUAD
        @test Ferrite.nodes_to_vtkorder(cell) == [1, 2, 3, 4]
        # 3D: bottom face followed by top face, as for the VTK wedge and hexahedron
        cell = InterfaceCell(Triangle((1, 2, 3)), Triangle((4, 5, 6)))
        @test Ferrite.cell_to_vtkcell(typeof(cell)) == Ferrite.VTKCellTypes.VTK_WEDGE
        @test Ferrite.nodes_to_vtkorder(cell) == [1, 2, 3, 4, 5, 6]
        cell = InterfaceCell(Quadrilateral((1, 2, 3, 4)), Quadrilateral((5, 6, 7, 8)))
        @test Ferrite.cell_to_vtkcell(typeof(cell)) == Ferrite.VTKCellTypes.VTK_HEXAHEDRON
        @test Ferrite.nodes_to_vtkorder(cell) == [1, 2, 3, 4, 5, 6, 7, 8]
    end

    @testset "Construction from node tuple" begin
        for cell in (
                InterfaceCell(Line((1, 2)), Line((4, 3))),
                InterfaceCell(QuadraticLine((1, 2, 5)), QuadraticLine((4, 3, 6))),
                InterfaceCell(Triangle((1, 2, 3)), Triangle((4, 5, 6))),
                InterfaceCell(QuadraticTriangle((1, 2, 3, 7, 8, 9)), QuadraticTriangle((4, 5, 6, 10, 11, 12))),
                InterfaceCell(Quadrilateral((1, 2, 3, 4)), Quadrilateral((5, 6, 7, 8))),
            )
            @test typeof(cell)(cell.nodes) == cell
        end
    end

    @testset "Write mixed grid" begin
        mktempdir() do tmp
            for grid in [
                    generate_grid(Triangle, (4, 4)),
                    generate_grid(Quadrilateral, (4, 4)),
                    generate_grid(Tetrahedron, (4, 4, 4)),
                    generate_grid(Hexahedron, (4, 4, 4)),
                ]
                addcellset!(grid, "A", x -> norm(x, Inf) ≤ 0.5)
                addcellset!(grid, "B", setdiff(OrderedSet(1:getncells(grid)), getcellset(grid, "A")))
                grid2 = insert_interfaces(grid, ["A", "B"])
                # The VTK cell of a zero-thickness interface must not be a bow-tie: nodes
                # crossing the interface (2 -> 3 and 4 -> 1 for the quad) must coincide.
                if Ferrite.getspatialdim(grid) == 2
                    for i in getcellset(grid2, "interfaces")
                        x = map(n -> get_node_coordinate(grid2, n), Ferrite.nodes_to_vtkorder(getcells(grid2, i)))
                        @test x[2] ≈ x[3] && x[4] ≈ x[1] && !(x[1] ≈ x[2])
                    end
                end
                dh = DofHandler(grid2)
                set_bulk = union(getcellset(grid2, "A"), getcellset(grid2, "B"))
                add!(SubDofHandler(dh, set_bulk), :u, Lagrange{Ferrite.getrefshape(grid.cells[1]), 1}())
                set_interface = getcellset(grid2, "interfaces")
                add!(SubDofHandler(dh, set_interface), :u, InterfaceCellInterpolation(Lagrange{FerriteInterfaceElements.getinterfaceshape(grid.cells[1]), 1}()))
                close!(dh)
                u = rand(ndofs(dh))
                for discontinuous in (false, true)
                    file = joinpath(tmp, "output_$discontinuous")
                    VTKGridFile(file, grid2; write_discontinuous = discontinuous) do vtk
                        write_solution(vtk, dh, u)
                        Ferrite.write_cellset(vtk, grid2)
                    end
                    @test isfile(file * ".vtu")
                end
            end
        end
    end
end
