function test_grid_data(old::Grid, new::Grid)
        # Check if all InterfaceCell's contain new nodes on at least one side and that node pairs have the same coordinates
    n = length(old.nodes)
    nodepairs = Set{Tuple{Int,Int}}()
    interfacecellids = Set{Int}()
    for (i, cell) in enumerate(new.cells)
        if cell isa InterfaceCell
            push!(interfacecellids, i)
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
        # Reconstruct cellsets based on interfaces
    nodesets = [Set{Int}(), Set{Int}(), Set{Int}(), Set{Int}()]
    cell = new.cells[pop!(interfacecellids)]
    for (n, m) in zip(cell.here.nodes, cell.there.nodes)
        push!(nodesets[1], n)
        push!(nodesets[2], m)
    end
            # Find nodes on same side of interfaces
    while ! isempty(interfacecellids)
        cell = new.cells[pop!(interfacecellids)]
        for nodes in (cell.here.nodes, cell.there.nodes)
                # Search for a set with node overlap
            foundconnection = false
            for set in nodesets
                if any( [n in set for n in nodes] )
                    foundconnection = true
                    for n in nodes
                        push!(set, n)
                    end
                    break
                end
            end
                # If no node overlap, add nodes to empty set
            if ! foundconnection
                for set in nodesets
                    if isempty(set)
                        for n in nodes
                            push!(set, n)
                        end
                        break
                    end
                end
            end
                # Merge sets with overlap
            for i in 1:length(nodesets)
                for j in i+1:length(nodesets)
                    if ! isempty( nodesets[i] ∩ nodesets[j] )
                        union!(nodesets[i], nodesets[j])
                        empty!(nodesets[j])
                    end
                end
            end
        end
    end
    filter!(s -> ! isempty(s), nodesets) # One set should be empty
        # Collect cells according to node sets
    cellsets = [OrderedSet{Int}(), OrderedSet{Int}(), OrderedSet{Int}()]
    for (cellid, cell) in enumerate(new.cells)
        if ! (cell isa InterfaceCell)
            for (i, set) in enumerate(nodesets)
                if any([n in set for n in cell.nodes])
                    push!(cellsets[i], cellid)
                    break
                end
            end
        end
    end
        # Test reconstructed cellsets
    for set in cellsets
        @test any( [sort(set) == sort(oldset) for oldset in values(old.cellsets)] )
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

@testset "Inserting interfaces in 2D for a mixed grid" begin
    # 7 ___ 8 ___ 9              7 ___ 8__14___ 9
    # |     |     |              |     |  |     |
    # | c4  | c5  |              |  c4 |c8| c5  |
    # |     |     |              |     |  |     |
    # 4 ___ 5 ___ 6     --->     11___10__13___12 
    # |\ c2 |     |              | c6   \/ c7 |
    # |  \  | c3  |              4 ___  5  ___ 6
    # | c1 \|     |              |\ c2  |      |
    # 1 ___ 2 ___ 3              |  \   |  c3  |
    #                            | c1 \ |      |
    #                            1 ____ 2 ____ 3
    #
    nodes = Node.([ Vec((-1.0, -1.0)), Vec((0.0, -1.0)), Vec((1.0, -1.0)),
                    Vec((-1.0,  0.0)), Vec((0.0,  0.0)), Vec((1.0,  0.0)),
                    Vec((-1.0,  1.0)), Vec((0.0,  1.0)), Vec((1.0,  1.0))])
    cells = [Triangle((1,2,4)), Triangle((2,5,4)), Quadrilateral((2,3,6,5)), Quadrilateral((4,5,8,7)), Quadrilateral((5,6,9,8))]
    grid = Grid(cells, nodes)
    addcellset!(grid, "bottom", OrderedSet((1,2,3)))
    addcellset!(grid, "topleft", OrderedSet((4,)))
    addcellset!(grid, "topright", OrderedSet((5,)))

    domain_names = ["bottom", "topleft", "topright"]
    new_grid = insert_interfaces(grid, domain_names)
    test_grid_data(grid, new_grid)
end

@testset "Inserting interfaces in 3D for a mixed grid" begin
    nodes = Node.([ Vec((-1.0, -1.0, -1.0)), Vec(( 0.0, -1.0, -1.0)), Vec(( 0.0, 0.0, -1.0)), Vec((-1.0,  0.0, -1.0)), # bottom of pyramids 1+2
                    Vec((-1.0, -1.0,  0.0)), # tip of pyramid 1
                    Vec((-1.0, -1.0, -2.0)), # tip of pyramid 2
                    Vec((-1.0, -2.0, -1.0)), Vec(( 0.0, -2.0, -1.0)), # bottom of pyramid 3
                    Vec((-1.0, -1.0, -2.0)), Vec((-1.0, -2.0, -2.0)), # wedge 1
                    ])
    cells = [Pyramid((1,2,3,4,5)), Pyramid((1,2,3,4,6)), Pyramid((1,2,8,7,5)), Wedge((1,2,9,8,7,10))] 
    grid = Grid(cells, nodes)
    addcellset!(grid, "p 1", OrderedSet((1,)))
    addcellset!(grid, "p 2", OrderedSet((2,)))
    addcellset!(grid, "p 3", OrderedSet((3,)))
    addcellset!(grid, "w 1", OrderedSet((4,)))
    
    domain_names = ["p 1", "p 2", "p 3", "w 1"]
    new_grid = insert_interfaces(grid, domain_names)
    @test length(new_grid.cells) == 7
    @test length(new_grid.nodes) == 21
end