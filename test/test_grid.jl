function test_grid_data(old::Grid, new::Grid)
        # Check if all InterfaceCell's contain new nodes on at least one side and that node pairs have the same coordinates
    n = length(old.nodes)
    nodepairs = Set{Tuple{Int,Int}}()
    for cell in new.cells
        if cell isa InterfaceCell
            h, t = cell.here.nodes, cell.there.nodes
            @test (any( h .≤ n ) && all( t .> n )) || (all( h .> n ) && any( t .≤ n ))
            for (n, m) in zip(h, t)
                push!(nodepairs, (n,m))
            end
        end
    end
    for (n, m) in nodepairs
        @test new.nodes[n] == new.nodes[m]
    end
    return nothing
end

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
    test_grid_data(grid, new_grid)
end

@testset "Inserting interfaces in 3D" begin
    grid = generate_grid(Hexahedron, (2,2,1))
    addcellset!(grid, "bottomleft", OrderedSet((1,)))
    addcellset!(grid, "topleft", OrderedSet((3,)))
    addcellset!(grid, "right", OrderedSet((2,4)))
    
    domain_names = ["bottomleft", "topleft", "right"]
    new_grid = insert_interfaces(grid, domain_names)
    test_grid_data(grid, new_grid)

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
