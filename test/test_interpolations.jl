@testset "InterfaceCellInterpolation" begin
    for (iphere, ipthere, shape) in (
            (Lagrange{RefTriangle, 1}(), Lagrange{RefTriangle, 1}(), RefPrism), 
            (Lagrange{RefTriangle, 2}(), Lagrange{RefTriangle, 2}(), RefPrism),
            (Lagrange{RefTriangle, 1}(), Lagrange{RefTriangle, 2}(), RefPrism),
            (Lagrange{RefQuadrilateral, 1}(), Lagrange{RefQuadrilateral, 1}(), RefHexahedron), )
        IPhere  = typeof(iphere)
        IPthere = typeof(ipthere)

        @test InterfaceCellInterpolation(iphere, ipthere) isa InterfaceCellInterpolation{shape, IPhere, IPthere}
        @test InterfaceCellInterpolation(iphere) isa InterfaceCellInterpolation{shape, IPhere, IPhere}
        ip = InterfaceCellInterpolation(iphere, ipthere)
        @test Ferrite.nvertices(ip) == Ferrite.nvertices(iphere) + Ferrite.nvertices(ipthere)
        @test all(Ferrite.vertexdof_indices(ip) .== collect( (v,) for v in 1:Ferrite.nvertices(ip) ))
    end
    here, there = Lagrange{RefTriangle, 1}(), Lagrange{RefTriangle, 2}()
    interface = InterfaceCellInterpolation(here, there)

    @test getnbasefunctions(interface) == 9

    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(interface, :here, 10)
    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(interface, :there, 10)
    @test_throws ArgumentError FerriteInterfaceElements.get_interface_index(interface, :nowhere, 10)

    @test Ferrite.facedof_indices(interface) == ((1,2,3), (4,5,6,7,8,9))

    @test Ferrite.edgedof_indices(interface) == ((1,2),(2,3),(3,1),(4,5,7),(5,6,8),(6,4,9))

    testcelltype = InterfaceCell{RefQuadrilateral, Line, QuadraticLine}
    expectedtype = InterfaceCellInterpolation{RefQuadrilateral, Lagrange{RefLine,1,Nothing}, Lagrange{RefLine,2,Nothing}}
    @test Ferrite.default_interpolation(testcelltype) isa expectedtype
    @test Ferrite.default_geometric_interpolation(Ferrite.default_interpolation(testcelltype)) isa VectorizedInterpolation{2, RefQuadrilateral, <:Any, expectedtype}
end
