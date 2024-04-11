@testset "Inserting interfaces in 2D" begin
    grid = generate_grid(Triangle, (2,2), Vec((-1.0,-1.0)), Vec((1.0,1.0)))
    
    #  7----8----9
    #  | \6 | \8 |
    #  | 5\ | 7\ |
    #  4----5----6
    #  | \2 | \4 |
    #  | 1\ | 3\ |
    #  1----2----3

    addcellset!(grid, "bottom",  x -> x[2] ≤ 0)
    addcellset!(grid, "top",     x -> x[2] ≥ 0)
    addcellset!(grid, "top left",  x -> x[1] ≤ 0 && x[2] ≥ 0)
    addcellset!(grid, "top right", x -> x[1] ≥ 0 && x[2] ≥ 0)

    newcells = create_interface_cells!(grid, "bottom", "top")
    @test newcells[1].nodes == (12,10,6,5) # bottom: 4 -> 11, 5 -> 10, 6 -> 12
    @test newcells[2].nodes == (10,11,5,4)
    for (i, n) in enumerate(((1,2,11), (2,10,11), (2,3,10), (3,12,10), (4,5,7), (5,8,7), (5,6,8), (6,9,8)))
        @test grid.cells[i].nodes == n
    end
    grid = Grid(vcat(grid.cells, newcells), grid.nodes; cellsets=grid.cellsets)

    newcells = create_interface_cells!(grid, "top left", "top right")
    @test newcells[1].nodes == (13,14,5,8) # top left: 5 -> 13, 8 -> 14
    for (i, n) in enumerate(((1,2,11), (2,10,11), (2,3,10), (3,12,10), (4,13,7), (13,14,7), (5,6,8), (6,9,8), (12,10,6,5), (10,11,13,4)))
        @test grid.cells[i].nodes == n
    end
end
