@testset "InterfaceCellInterpolation" begin
    for (baseshape, shape) in ( (RefLine, RefQuadrilateral), (RefTriangle, RefPrism), (RefQuadrilateral, RefHexahedron) )
        for order in (1,2)
            IP   = Lagrange{baseshape, order}
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
    @test Ferrite.facetdof_interior_indices(ip) == ((17,), (18,))
    @test Ferrite.edgedof_interior_indices(ip) == ((9,),(10,),(11,),(12,),(13,),(14,),(15,),(16,))

    @test get_side_and_baseindex(ip, 5) == (:there, 1)
    @test_throws ArgumentError get_side_and_baseindex(ip, 19)

    testcelltype = InterfaceCell{RefQuadrilateral, Line, 4}
    expectedtype = InterfaceCellInterpolation{RefQuadrilateral, 1, Lagrange{RefLine,1}}
    @test Ferrite.default_geometric_interpolation(testcelltype) isa expectedtype
    @test Ferrite.default_geometric_interpolation(Ferrite.default_geometric_interpolation(testcelltype)) isa VectorizedInterpolation{2, RefQuadrilateral, <:Any, expectedtype}

    @test_throws AssertionError FerriteInterfaceElements.get_interface_index(ip, :here, 0)
    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(ip, :here, 100)
    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(ip, :there, 100)
    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(ip, :test, 1)

    @test length(Ferrite.reference_coordinates(ip)) == getnbasefunctions(ip)
end

@testset "Higher-order interface DOF numbering" begin
    for base in (Lagrange{RefTriangle,3}(), Lagrange{RefTriangle,4}(), Lagrange{RefQuadrilateral,3}())
        ip = InterfaceCellInterpolation(base)
        vertices = collect(Iterators.flatten(Ferrite.vertexdof_indices(ip)))
        edges = collect(Iterators.flatten(Ferrite.edgedof_interior_indices(ip)))
        faces = collect(Iterators.flatten(Ferrite.facedof_interior_indices(ip)))
        @test vcat(vertices, edges, faces) == collect(1:getnbasefunctions(ip))
        for entitydofs in (Ferrite.edgedof_interior_indices, Ferrite.facedof_interior_indices)
            base_entities = entitydofs(base)
            interface_entities = entitydofs(ip)
            n = length(base_entities)
            for j in 1:n
                @test get_side_and_baseindex.(Ref(ip), interface_entities[j]) == map(i -> (:here, i), base_entities[j])
                @test get_side_and_baseindex.(Ref(ip), interface_entities[n+j]) == map(i -> (:there, i), base_entities[j])
            end
        end
    end
end
