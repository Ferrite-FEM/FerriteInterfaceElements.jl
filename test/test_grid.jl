@testset "Inserting interfaces in 2D" begin
    # 7 ___ 8 ___ 9              7 ___ 8__14___ 9
    # |\ c6 |\ c8 |              |\ c6 | c|\ c8 |
    # |  \  |  \  |              |  \  |11|  \  |
    # | c5 \| c7 \|              | c5 \|  | c7 \|
    # 4 ___ 5 ___ 6     --->     13___12__11___10 
    # |\ c2 |\ c4 |              | c10  \/  c9 |
    # |  \  |  \  |              4 ___  5  ___ 6
    # | c1 \| c3 \|              |\ c2  | \ c4 |
    # 1 ___ 2 ___ 3              |  \   |   \  |
    #                            | c1 \ |  c3 \|
    #                            1 ____ 2 ____ 3
    #
    grid = generate_grid(Triangle, (2,2))
    addcellset!(grid, "bottom", Set((1,2,3,4)))
    addcellset!(grid, "topleft", Set((5,6)))
    addcellset!(grid, "topright", Set((7,8)))

    domain_names = ["bottom", "topleft", "topright"]
    new_grid = insert_interfaces(grid, domain_names)

    for (nodeid, duplicate_nodeid) in [(4,13), (5,12), (5,11), (11,12), (6,10), (8,14)]
        @test new_grid.nodes[nodeid] == new_grid.nodes[duplicate_nodeid]
    end

    @test new_grid.cells == [
        Triangle((1, 2, 4)),
        Triangle((2, 5, 4)),
        Triangle((2, 3, 5)),
        Triangle((3, 6, 5)),
        Triangle((13, 12, 7)),
        Triangle((12, 8, 7)),
        Triangle((11, 10, 14)),
        Triangle((10, 9, 14)),
        InterfaceCell(Line((6, 5)), Line((10, 11))),
        InterfaceCell(Line((5, 4)), Line((12, 13))),
        InterfaceCell(Line((12, 8)), Line((11, 14)))
        ]
end

@testset "Inserting interfaces in 3D" begin
    grid = generate_grid(Hexahedron, (2,2,1))
    addcellset!(grid, "bottomleft", Set((1,)))
    addcellset!(grid, "topleft", Set((3,)))
    addcellset!(grid, "right", Set((2,4)))
    
    domain_names = ["bottomleft", "topleft", "right"]
    new_grid = insert_interfaces(grid, domain_names)

    for (nodeid, duplicate_nodeid) in [(13, 25), (14,26), (14,21), (21,26), (11,22), (17,27), (4,24), (5,23), (5,20), (20,23), (8, 28), (2,19)]
        @test new_grid.nodes[nodeid] == new_grid.nodes[duplicate_nodeid]
    end

    @test new_grid.cells ==  [
        Hexahedron((1, 2, 5, 4, 10, 11, 14, 13)),
        Hexahedron((19, 3, 6, 20, 22, 12, 15, 21)),
        Hexahedron((24, 23, 28, 7, 25, 26, 27, 16)),
        Hexahedron((20, 6, 9, 8, 21, 15, 18, 17)),
        InterfaceCell(Quadrilateral((2, 5, 14, 11)), Quadrilateral((19, 20, 21, 22))),
        InterfaceCell(Quadrilateral((5, 4, 13, 14)), Quadrilateral((23, 24, 25, 26))),
        InterfaceCell(Quadrilateral((20, 21, 17, 8)), Quadrilateral((23, 26, 27, 28)))
        ]

    # More complex example with rough test
    grid = generate_grid(Tetrahedron, (10,10,10), Vec((-0.5,-0.5,-0.5)), Vec((0.5,0.5,0.5)))
    addcellset!(grid, "1", x -> true)
    set1 = getcellset(grid, "1")
    addcellset!(grid, "2", x -> (-0.3 ≤ x[1] ≤ 0.3) &&
                                (-0.3 ≤ x[2] ≤ 0.3) &&
                                (-0.3 ≤ x[3] ≤ 0.3))
    set2 = getcellset(grid, "2")
    addcellset!(grid, "3", x -> ((-0.1 ≤ x[1] ≤ 0.1) && (-0.1 ≤ x[2] ≤ 0.1)) ||
                                ((-0.1 ≤ x[1] ≤ 0.1) && (-0.1 ≤ x[3] ≤ 0.1)) ||
                                ((-0.1 ≤ x[2] ≤ 0.1) && (-0.1 ≤ x[3] ≤ 0.1)))
    set3 = getcellset(grid, "3")
    setdiff!(set1, set2, set3)
    setdiff!(set3, set2)

    newgrid = insert_interfaces(grid, ["1", "2", "3"])
    @test length(getcells(newgrid, "1-2-interface")) == 384
    @test length(getcells(newgrid, "1-3-interface")) == 192
    @test length(getcells(newgrid, "2-3-interface")) == 48
end
