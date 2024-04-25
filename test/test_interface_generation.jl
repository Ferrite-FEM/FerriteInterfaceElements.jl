##################################
# test stuff
#
# 7 ___ 8 ___ 9
# |\ c6 |\ c8 |
# |  \  |  \  |
# | c5 \| c7 \|
# 4 ___ 5 ___ 6
# |\ c2 |\ c4 |
# |  \  |  \  |
# | c1 \| c3 \|
# 1 ___ 2 ___ 3
#
# 7 ___ 8__14___ 9
# |\ c6 | c|\ c8 |
# |  \  |11|  \  |
# | c5 \|  | c7 \|
# 13___12__11___10
# | c10  \/  c9 |
# 4 ___  5  ___ 6
# |\ c2  | \ c4 |
# |  \   |   \  |
# | c1 \ |  c3 \|
# 1 ____ 2 ____ 3



grid = generate_grid(Triangle, (2,2))
addcellset!(grid, "bottom", Set((1,2,3,4)))
addcellset!(grid, "topleft", Set((5,6)))
addcellset!(grid, "topright", Set((7,8)))

topology = ExclusiveTopology(grid)

domain_names = ["bottom", "topleft", "topright"]

new_grid = FerriteCohesiveZones.insert_interfaces(grid, domain_names; topology)

for (nodeid, duplicate_nodeid) in [(4,13), (5,12), (5,11), (11,12), (6,10), (8,14)]
    @test new_grid.nodes[nodeid] == new_grid.nodes[duplicate_nodeid]
end

new_grid.cells == [ Triangle((1, 2, 4)),
                    Triangle((2, 5, 4)),
                    Triangle((2, 3, 5)),
                    Triangle((3, 6, 5)),
                    Triangle((13, 12, 7)),
                    Triangle((12, 8, 7)),
                    Triangle((11, 10, 14)),
                    Triangle((10, 9, 14)),
                    CohesiveQuadrilateral((6, 5, 10, 11)),
                    CohesiveQuadrilateral((5, 4, 12, 13)),
                    CohesiveQuadrilateral((12, 8, 11, 14)),]

# 3D 
grid = generate_grid(Hexahedron, (2,2,1))
addcellset!(grid, "bottomleft", Set((1,)))
addcellset!(grid, "topleft", Set((3,)))
addcellset!(grid, "right", Set((2,4)))
domain_names = ["bottomleft", "topleft", "right"]

new_grid = FerriteCohesiveZones.insert_interfaces(grid, domain_names)

@test grid.cells[1] == new_grid.cells[1]
@test new_grid.cells ==  [   Hexahedron((1, 2, 5, 4, 10, 11, 14, 13)),
                         Hexahedron((19, 3, 6, 20, 22, 12, 15, 21)),
                         Hexahedron((24, 23, 28, 7, 25, 26, 27, 16)),
                         Hexahedron((20, 6, 9, 8, 21, 15, 18, 17)),
                         FerriteCohesiveZones.CohesiveHexahedron((2, 5, 14, 11, 19, 20, 21, 22)),
                         FerriteCohesiveZones.CohesiveHexahedron((5, 4, 13, 14, 23, 24, 25, 26)),
                         FerriteCohesiveZones.CohesiveHexahedron((20, 21, 17, 8, 23, 26, 27, 28)),]

for (nodeid, duplicate_nodeid) in [(13, 25), (14,26), (14,21), (21,26), (11,22), (17,27), (4,24), (5,23), (5,20), (20,23), (8, 28), (2,19)]
    @test new_grid.nodes[nodeid] == new_grid.nodes[duplicate_nodeid]
end
