@testset "InterfaceCell" begin
    for (here, there, shape) in (
            (         Line((1,2)),            Line((3,4)),   RefQuadrilateral), 
            (QuadraticLine((1,2,5)), QuadraticLine((3,4,6)), RefQuadrilateral),
            (         Triangle((1,2,3)),                Triangle((4,5,6)),          RefPrism), 
            (QuadraticTriangle((1,2,3,7,8,9)), QuadraticTriangle((4,5,6,10,11,12)), RefPrism),
            (         Quadrilateral((1,2,3,4)),                        Quadrilateral((5,6,7,8)),                RefHexahedron), 
            (QuadraticQuadrilateral((1,2,3,4,9,10,11,12,17)), QuadraticQuadrilateral((5,6,7,8,13,14,15,16,18)), RefHexahedron)
            ) # The nodes have been chosen as the numbers representing their expected ordering
        C  = typeof(here)

        @test InterfaceCell{shape, C}(here, there) isa InterfaceCell{shape, C}
        @test InterfaceCell(here, there) isa InterfaceCell{shape, C}
        cell = InterfaceCell(here, there)

        @test Ferrite.nvertices(cell) == Ferrite.nvertices(here) + Ferrite.nvertices(there)
        @test Ferrite.nfaces(cell) == 2
        @test Ferrite.nnodes(cell) == Ferrite.nnodes(here) + Ferrite.nnodes(there)
        
        @test Ferrite.get_node_ids(cell) == ntuple(i -> i, Ferrite.nnodes(cell))
        @test Ferrite.vertices(cell) == (Ferrite.vertices(here)..., Ferrite.vertices(there)...)
        @test Ferrite.faces(cell) == (Ferrite.vertices(here), Ferrite.vertices(there))
    end
end
