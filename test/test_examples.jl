function prepare_interface_test_grid_2D(order::Integer)
    nodes =[Node((0.0,0.0)), Node((1.0,0.0)), Node((1.0,1.0)), Node((0.0,1.0)),
            Node((0.5,0.0)), Node((1.0,0.5)), Node((0.5,1.0)), Node((0.0,0.5)),
            Node((0.5,0.5)),
            # <- left quadrilateral | right quadrilateral ->
            Node((1.0,0.0)), Node((2.0,0.0)), Node((2.0,1.0)), Node((1.0,1.0)),
            Node((1.5,0.0)), Node((2.0,0.5)), Node((1.5,1.0)), Node((1.0,0.5)),
            Node((1.5,0.5))
            ]
    if order == 1
        cells = [ 
                Quadrilateral( Tuple(1:4) ),
                InterfaceCell( Line( (2,3) ), Line( (10,13) ) ),
                Quadrilateral( Tuple(10:13) )
                ]
    elseif order == 2
        cells = [ 
                QuadraticQuadrilateral( Tuple(1:9) ),
                InterfaceCell(QuadraticLine( (2,3,6) ), QuadraticLine( (10,13,17) )),
                QuadraticQuadrilateral( Tuple(10:18) )
                ]
    end
    grid = Grid(cells, nodes)
    addfaceset!(grid, "∂Ωₗ", x -> x[1] ≈ 0)
    addfaceset!(grid, "∂Ωᵣ", x -> x[1] ≈ 2)
    return grid
end

function prepare_interface_test_grid_3D(order::Integer)
    nodes =[Node((0.0,0.0,0.0)), Node((1.0,0.0,0.0)), Node((1.0,1.0,0.0)), Node((0.0,1.0,0.0)),
            Node((0.0,0.0,1.0)), Node((1.0,0.0,1.0)), Node((1.0,1.0,1.0)), Node((0.0,1.0,1.0)),
            Node((0.5,0.0,0.0)), Node((1.0,0.5,0.0)), Node((0.5,1.0,0.0)), Node((0.0,0.5,0.0)), 
            Node((0.5,0.0,1.0)), Node((1.0,0.5,1.0)), Node((0.5,1.0,1.0)), Node((0.0,0.5,1.0)), 
            Node((0.0,0.0,0.5)), Node((1.0,0.0,0.5)), Node((1.0,1.0,0.5)), Node((0.0,1.0,0.5)), 
            Node((0.5,0.5,0.0)), Node((0.5,0.0,0.5)), Node((1.0,0.5,0.5)), Node((0.5,1.0,0.5)), Node((0.0,0.5,0.5)), Node((0.5,0.5,1.0)),
            Node((0.5,0.5,0.5)),
            # <- left hexahdron | right hexahedron ->
            Node((1.0,0.0,0.0)), Node((2.0,0.0,0.0)), Node((2.0,1.0,0.0)), Node((1.0,1.0,0.0)),
            Node((1.0,0.0,1.0)), Node((2.0,0.0,1.0)), Node((2.0,1.0,1.0)), Node((1.0,1.0,1.0)),
            Node((1.5,0.0,0.0)), Node((2.0,0.5,0.0)), Node((1.5,1.0,0.0)), Node((1.0,0.5,0.0)), 
            Node((1.5,0.0,1.0)), Node((2.0,0.5,1.0)), Node((1.5,1.0,1.0)), Node((1.0,0.5,1.0)), 
            Node((1.0,0.0,0.5)), Node((2.0,0.0,0.5)), Node((2.0,1.0,0.5)), Node((1.0,1.0,0.5)), 
            Node((1.5,0.5,0.0)), Node((1.5,0.0,0.5)), Node((2.0,0.5,0.5)), Node((1.5,1.0,0.5)), Node((1.0,0.5,0.5)), Node((1.5,0.5,1.0)),
            Node((1.5,0.5,0.5)),
            ]
    if order == 1
        cells = [ 
                Hexahedron( Tuple(1:8) ),
                InterfaceCell( Quadrilateral( (2,3,7,6) ), Quadrilateral( (28,31,35,32) ) ),
                Hexahedron( Tuple(28:35) )
                ]
    elseif order == 2
        cells = [ 
                QuadraticHexahedron( Tuple(1:27) ),
                InterfaceCell(QuadraticQuadrilateral( (2,3,7,6,10,19,14,18,23) ), QuadraticQuadrilateral( (28,31,35,32,39,47,43,44,52) )),
                QuadraticHexahedron( Tuple(28:54) )
                ]
    end
    grid = Grid(cells, nodes)
    addfaceset!(grid, "∂Ωₗ", x -> x[1] ≈ 0)
    addfaceset!(grid, "∂Ωᵣ", x -> x[1] ≈ 2)
    return grid
end

function prepare_interface_test_scalar_interpolation_2D(order::Integer)
    ip = (  left = Lagrange{RefQuadrilateral, order}(),
            interface = InterfaceCellInterpolation(Lagrange{RefLine, order}()),
            right = Lagrange{RefQuadrilateral, order}())

    cv = (  left = CellValues(QuadratureRule{RefQuadrilateral}(4), ip.left, ip.left),
            interface = InterfaceCellValues(QuadratureRule{RefLine}(4), ip.interface),
            right = CellValues(QuadratureRule{RefQuadrilateral}(4), ip.right, ip.right))
    return ip, cv
end

function prepare_interface_test_scalar_interpolation_3D(order::Integer)
    ip = (  left = Lagrange{RefHexahedron, order}(),
            interface = InterfaceCellInterpolation(Lagrange{RefQuadrilateral, order}()),
            right = Lagrange{RefHexahedron, order}())

    cv = (  left = CellValues(QuadratureRule{RefHexahedron}(4), ip.left, ip.left),
            interface = InterfaceCellValues(QuadratureRule{RefQuadrilateral}(4), ip.interface),
            right = CellValues(QuadratureRule{RefHexahedron}(4), ip.right, ip.right))
    return ip, cv
end

function prepare_interface_test_vector_interpolation_2D(order::Integer)
    ip = (  left = Lagrange{RefQuadrilateral, order}()^2,
            interface = InterfaceCellInterpolation(Lagrange{RefLine, order}())^2,
            right = Lagrange{RefQuadrilateral, order}()^2)

    cv = (  left = CellValues(QuadratureRule{RefQuadrilateral}(4), ip.left, ip.left),
            interface = InterfaceCellValues(QuadratureRule{RefLine}(4), ip.interface),
            right = CellValues(QuadratureRule{RefQuadrilateral}(4), ip.right, ip.right))
    return ip, cv
end

function prepare_interface_test_vector_interpolation_3D(order::Integer)
    ip = (  left = Lagrange{RefHexahedron, order}()^3,
            interface = InterfaceCellInterpolation(Lagrange{RefQuadrilateral, order}())^3,
            right = Lagrange{RefHexahedron, order}()^3)

    cv = (  left = CellValues(QuadratureRule{RefHexahedron}(4), ip.left, ip.left),
            interface = InterfaceCellValues(QuadratureRule{RefQuadrilateral}(4), ip.interface),
            right = CellValues(QuadratureRule{RefHexahedron}(4), ip.right, ip.right))
    return ip, cv
end

@testset "Example problems" begin
    @testset "Diffusion" begin
        include("test_example_diffusion.jl")
    end
    @testset "Cohesion" begin
        include("test_example_cohesion.jl")
    end
end 