@testset "InterfaceCellInterpolation" begin
    for (baseshape, shape) in ( (RefLine, RefQuadrilateral), (RefTriangle, RefPrism), (RefQuadrilateral, RefHexahedron) )
        for order in (1,2)
            IP   = Lagrange{baseshape, order, Nothing}
            base = Lagrange{baseshape, order}()
            @test InterfaceCellInterpolation(base) isa InterfaceCellInterpolation{shape, order, IP}
            ip = InterfaceCellInterpolation(base)
            @test Ferrite.nvertices(ip) == 2*Ferrite.nvertices(base)
            @test Ferrite.vertexdof_indices(ip) == Tuple( (v,) for v in 1:Ferrite.nvertices(ip) )
            @test Ferrite.getorder(ip) == order
        end
    end
    base = Lagrange{RefQuadrilateral, 2}()
    ip = InterfaceCellInterpolation(base)

    @test getnbasefunctions(ip) == 18

    @test Ferrite.vertexdof_indices(ip) == ((1,),(2,),(3,),(4,),(5,),(6,),(7,),(8,))
    @test Ferrite.facedof_interior_indices(ip) == ((17,), (18,))
    @test Ferrite.edgedof_interior_indices(ip) == ((9,),(10,),(11,),(12,),(13,),(14,),(15,),(16,))

    @test get_side_and_baseindex(ip, 5) == (:there, 1)
    @test_throws ArgumentError get_side_and_baseindex(ip, 19)

    testcelltype = InterfaceCell{RefQuadrilateral, Line, 4}
    expectedtype = InterfaceCellInterpolation{RefQuadrilateral, 1, Lagrange{RefLine,1,Nothing}}
    @test Ferrite.default_interpolation(testcelltype) isa expectedtype
    @test Ferrite.default_geometric_interpolation(Ferrite.default_interpolation(testcelltype)) isa VectorizedInterpolation{2, RefQuadrilateral, <:Any, expectedtype}
end
