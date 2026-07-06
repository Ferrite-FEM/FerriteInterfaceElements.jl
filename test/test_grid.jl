@testset "Inserting interfaces in 2D" begin
    # 7 ___ 8 ___ 9              7 ___ 8__14___ 9
    # |\ c6 |\ c8 |              |\ c6 | c|\ c8 |
    # |  \  |  \  |              |  \  |11|  \  |
    # | c5 \| c7 \|              | c5 \|  | c7 \|
    # 4 ___ 5 ___ 6     --->     11___10__13___12 
    # |\ c2 |\ c4 |              | c9   \/ c10 |
    # |  \  |  \  |              4 ___  5  ___ 6
    # | c1 \| c3 \|              |\ c2  | \ c4 |
    # 1 ___ 2 ___ 3              |  \   |   \  |
    #                            | c1 \ |  c3 \|
    #                            1 ____ 2 ____ 3
    #
    grid = generate_grid(Triangle, (2,2))
    addcellset!(grid, "bottom", OrderedSet((1,2,3,4)))
    addcellset!(grid, "topleft", OrderedSet((5,6)))
    addcellset!(grid, "topright", OrderedSet((7,8)))

    domain_names = ["bottom", "topleft", "topright"]
    new_grid = insert_interfaces(grid, domain_names)

    nodepairs = Set{Tuple{Int,Int}}()
    for cell in new_grid.cells
        if cell isa InterfaceCell
            for (n, m) in zip(cell.here.nodes, cell.there.nodes)
                push!(nodepairs, (n,m))
            end
        end
    end
    for (n, m) in nodepairs
        @test new_grid.nodes[n] == new_grid.nodes[m]
    end
    for cell in [Triangle((1, 2, 4)),
                 Triangle((2, 5, 4)),
                 Triangle((2, 3, 5)),
                 Triangle((3, 6, 5)),
                 Triangle((11, 10, 7)),
                 Triangle((10, 8, 7)),
                 Triangle((13, 12, 14)),
                 Triangle((12, 9, 14))]
        @test cell in new_grid.cells[1:8]
    end
    for cell in [InterfaceCell(Line((5, 4)), Line((10, 11))),
                 InterfaceCell(Line((6, 5)), Line((12, 13))),
                 InterfaceCell(Line((10, 8)), Line((13, 14)))]
        @test cell in new_grid.cells[9:11]
    end
end

@testset "Inserting interfaces in 3D" begin
    grid = generate_grid(Hexahedron, (2,2,1))
    addcellset!(grid, "bottomleft", OrderedSet((1,)))
    addcellset!(grid, "topleft", OrderedSet((3,)))
    addcellset!(grid, "right", OrderedSet((2,4)))
    
    domain_names = ["bottomleft", "topleft", "right"]
    new_grid = insert_interfaces(grid, domain_names)

    nodepairs = Set{Tuple{Int,Int}}()
    for cell in new_grid.cells
        if cell isa InterfaceCell
            for (n, m) in zip(cell.here.nodes, cell.there.nodes)
                push!(nodepairs, (n,m))
            end
        end
    end
    for (n, m) in nodepairs
        @test new_grid.nodes[n] == new_grid.nodes[m]
    end
    for cell in [Hexahedron((1, 2, 5, 4, 10, 11, 14, 13)),
                 Hexahedron((19, 3, 6, 20, 22, 12, 15, 21)),
                 Hexahedron((24, 23, 28, 7, 25, 26, 27, 16)),
                 Hexahedron((20, 6, 9, 8, 21, 15, 18, 17))]
        @test cell in new_grid.cells[1:4]
    end
    for cell in [InterfaceCell(Quadrilateral((2, 5, 14, 11)), Quadrilateral((19, 20, 21, 22))),
                 InterfaceCell(Quadrilateral((5, 4, 13, 14)), Quadrilateral((23, 24, 25, 26))),
                 InterfaceCell(Quadrilateral((20, 21, 17, 8)), Quadrilateral((23, 26, 27, 28)))]
        @test cell in new_grid.cells[5:7]
    end

    # More complex example with rough test
    grid = generate_grid(Tetrahedron, (10,10,10), Vec((-0.5,-0.5,-0.5)), Vec((0.5,0.5,0.5)))
    addcellset!(grid, "1", x -> true)
    set1 = getcellset(grid, "1")
    addcellset!(grid, "2", x -> (-0.31 < x[1] < 0.31) &&
                                (-0.31 < x[2] < 0.31) &&
                                (-0.31 < x[3] < 0.31))
    set2 = getcellset(grid, "2")
    addcellset!(grid, "3", x -> ((-0.11 < x[1] < 0.11) && (-0.11 < x[2] < 0.11)) ||
                                ((-0.11 < x[1] < 0.11) && (-0.11 < x[3] < 0.11)) ||
                                ((-0.11 < x[2] < 0.11) && (-0.11 < x[3] < 0.11)))
    set3 = getcellset(grid, "3")
    setdiff!(set1, set2, set3)
    setdiff!(set3, set2)

    newgrid = insert_interfaces(grid, ["1", "2", "3"])
    @test length(getcells(newgrid, "1-2-interface")) == 384
    @test length(getcells(newgrid, "1-3-interface")) == 192
    @test length(getcells(newgrid, "2-3-interface")) == 48
end
