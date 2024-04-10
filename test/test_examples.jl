function prepare_interface_test_grid(order::Tuple{<:Integer,<:Integer})
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
    cells = Ferrite.AbstractCell[]
    
    order[1] == 1 ? push!(cells, Hexahedron( (collect(1:8)...,) )) : nothing
    order[1] == 2 ? push!(cells, QuadraticHexahedron( (collect(1:27)...,) )) : nothing
    
    order == (1,1) ? push!(cells, InterfaceCell( Quadrilateral( (2,3,7,6) ), Quadrilateral( (28,31,35,32) ) )) : nothing
    order == (2,2) ? push!(cells, InterfaceCell(QuadraticQuadrilateral( (2,3,7,6,10,19,14,18,23) ), QuadraticQuadrilateral( (28,31,35,32,39,47,43,44,52) ))) : nothing
    order == (1,2) ? push!(cells, InterfaceCell( Quadrilateral( (2,3,7,6) ), QuadraticQuadrilateral( (28,31,35,32,39,47,43,44,52) ) )) : nothing
    order == (2,1) ? push!(cells, InterfaceCell(QuadraticQuadrilateral( (2,3,7,6,10,19,14,18,23) ), Quadrilateral( (28,31,35,32) ) )) : nothing
                  
    order[2] == 1 ? push!(cells, Hexahedron( (collect(28:35)...,) )) : nothing
    order[2] == 2 ? push!(cells, QuadraticHexahedron( (collect(28:54)...,) )) : nothing
    
    grid = Grid(cells, nodes)
    
    addfaceset!(grid, "∂Ωₗ", x -> x[1] ≈ 0)
    addfaceset!(grid, "∂Ωᵣ", x -> x[1] ≈ 2)
return grid
end
    
@testset "Test-example for scalar field" begin
    @testset "ip-order=$order" for order in ((1,1), (2,2), (1,2))  
        function assemble_test_element!(Kₑ::Matrix, cv::CellValues)
            nbf = getnbasefunctions(cv)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                for i in 1:nbf
                    ∇δc = shape_gradient(cv, qp, i)
                    for j in 1:nbf
                    ∇c = shape_gradient(cv, qp, j)
                        Kₑ[i,j] += (∇δc ⋅ ∇c) * dΩ
                    end
                end
            end
        end
        
        function assemble_test_element!(Kₑ::Matrix, cv::InterfaceCellValues)
            nbf = getnbasefunctions(cv)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv.here, qp)
                for i in 1:nbf
                    side, bi = get_side_and_baseindex(cv, i)
                    side == :here ? jump_δc = -shape_value(cv.here, qp, bi) : jump_δc = shape_value(cv.there, qp, bi)
                    for j in 1:nbf
                        side, bi = get_side_and_baseindex(cv, j)
                        side == :here ? jump_c = -shape_value(cv.here, qp, bi) : jump_c = shape_value(cv.there, qp, bi)
                        Kₑ[i,j] += (jump_δc * jump_c) * dΩ
                    end
                end
            end
        end
    
        grid = prepare_interface_test_grid(order)
        cellid = (left=1, interface=2, right=3)
    
        ip = (  left = Lagrange{RefHexahedron, order[1]}(),
                interface = InterfaceCellInterpolation(Lagrange{RefQuadrilateral, order[1]}(), Lagrange{RefQuadrilateral, order[2]}()), 
                right = Lagrange{RefHexahedron, order[2]}())
    
        cv = (  left = CellValues(QuadratureRule{RefHexahedron}(4), ip.left, ip.left), # Default geometry interpolation resulting in exception for order=2
                interface = InterfaceCellValues(QuadratureRule{RefQuadrilateral}(4), ip.interface), #, InterfaceCellInterpolation(Lagrange{RefQuadrilateral, order[1]}(), Lagrange{RefQuadrilateral, order[2]}()))
                right = CellValues(QuadratureRule{RefHexahedron}(4), ip.right, ip.right))
    
        dh = DofHandler(grid)
        function prepare_sub_dof_handler(key::Symbol)
            s = SubDofHandler(dh, Set{Int}([cellid[key]]))
            if order[1] == order[2] || key == :left
                add!(s, :c, ip[key])
            else
                warntext = "Field :c uses a different interpolation order in another SubDofHandler."
                if key == :interface
                    @test_logs (:warn, warntext) add!(s, :c, ip[key])
                else
                    @test_logs (:warn, warntext) (:warn, warntext) add!(s, :c, ip[key])
                end
            end
            return s
        end
        sdh = ( left = prepare_sub_dof_handler(:left),
                interface = prepare_sub_dof_handler(:interface),
                right = prepare_sub_dof_handler(:right))
        close!(dh)
    
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:c, getfaceset(grid, "∂Ωₗ"), (x,t) -> 3.0))
        add!(ch, Dirichlet(:c, getfaceset(grid, "∂Ωᵣ"), (x,t) -> 0.0))
        close!(ch)
    
        @testset "Accessing shape values etc." begin
            icv = cv[:interface]
            for cc in CellIterator(dh, Set{Int}([cellid[:interface]]))
                reinit!(icv, cc)
            end
            nbf = getnbasefunctions(icv)
            qp = rand(1:getnquadpoints(icv))
            i = rand(1:nbf)
            side, bi = get_side_and_baseindex(icv, i)
            if side == :here
                @test shape_value(icv, qp, i, true) == shape_value(icv.here, qp, bi)
                @test shape_value(icv, qp, i, false) == 0
                @test shape_value_average(icv, qp, i) == 0.5*shape_value(icv.here, qp, bi)
                @test shape_gradient_average(icv, qp, i) == 0.5*shape_gradient(icv.here, qp, bi)
                @test shape_value_jump(icv, qp, i) == -shape_value(icv.here, qp, bi)
                @test shape_gradient_jump(icv, qp, i) == -shape_gradient(icv.here, qp, bi)
            else
                @test shape_value(icv, qp, i, true) == 0
                @test shape_value(icv, qp, i, false) == shape_value(icv.there, qp, bi)
                @test shape_value_average(icv, qp, i) == 0.5*shape_value(icv.there, qp, bi)
                @test shape_gradient_average(icv, qp, i) == 0.5*shape_gradient(icv.there, qp, bi)
                @test shape_value_jump(icv, qp, i) == shape_value(icv.there, qp, bi)
                @test shape_gradient_jump(icv, qp, i) == shape_gradient(icv.there, qp, bi)
            end
        end
    
        K = create_sparsity_pattern(dh, ch)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K,f)
    
        for key in (:left, :interface, :right)
            nbf = getnbasefunctions(cv[key])
            Kₑ = zeros(nbf, nbf)
            for cc in CellIterator(dh, Set{Int}([cellid[key]]))
                fill!(Kₑ, 0)
                reinit!(cv[key], cc)
                assemble_test_element!(Kₑ, cv[key])
                assemble!(assembler, celldofs(cc), Kₑ)
            end
        end
    
        apply!(K, f, ch)
        u = K \ f
    
        expectedvalue = 3.0
        for (cellid, leftface, rightface) in ((1,5,3), (2,1,2), (3,5,3))
            cell = getcells(grid, cellid)
            cellfaces = Ferrite.faces(cell)
            cellnodes = Ferrite.get_node_ids(cell)
            dofs = celldofs(dh, cellid)
            leftdofs = dofs[indexin([cellfaces[leftface]...], [cellnodes...])]
            leftvalues = u[leftdofs]
            rightdofs = dofs[indexin([cellfaces[rightface]...], [cellnodes...])]
            rightvalues = u[rightdofs]
            @test all(leftvalues .≈ expectedvalue)
            @test all(rightvalues .≈ expectedvalue-1)
            expectedvalue -= 1
        end
    end
end
    
@testset "Test-example for vector field" begin
    @testset "ip-order=$order" for order in ((1,1), (2,2), (1,2))
        δ(i,j) = i == j ? 1.0 : 0.0
        E = SymmetricTensor{4,3}( (i,j,k,l) -> 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) ) # Using E=1, ν=0
    
        function assemble_test_element!(Kₑ::Matrix, cv::CellValues)
            nbf = getnbasefunctions(cv)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv, qp)
                for i in 1:nbf
                    δε   = shape_symmetric_gradient(cv, qp, i)
                    for j in 1:nbf
                        ε = shape_symmetric_gradient(cv, qp, j)
                        Kₑ[i,j] += (δε ⊡ E ⊡ ε) * dΩ
                    end
                end
            end
        end
        
        function assemble_test_element!(Kₑ::Matrix, cv::InterfaceCellValues)
            nbf = getnbasefunctions(cv)
            for qp in 1:getnquadpoints(cv)
                dΩ = getdetJdV(cv.here, qp)
                for i in 1:nbf
                    s, bi = get_side_and_baseindex(cv, i)
                    s == :here ? jump_δu = -shape_value(cv.here, qp, bi) : jump_δu = shape_value(cv.there, qp, bi)
                    for j in 1:nbf
                        s, bi = get_side_and_baseindex(cv, j)
                        s == :here ? jump_u = -shape_value(cv.here, qp, bi) : jump_u = shape_value(cv.there, qp, bi)
                        Kₑ[i,j] += 1e10*(jump_δu ⋅ jump_u) * dΩ # "Gluing together the two bulk cells"
                    end
                end
            end
        end
    
        grid = prepare_interface_test_grid(order)
        cellid = (left=1, interface=2, right=3)
    
        ip = (  left = Lagrange{RefHexahedron, order[1]}()^3,
                interface = InterfaceCellInterpolation(Lagrange{RefQuadrilateral, order[1]}(), Lagrange{RefQuadrilateral, order[2]}())^3,
                right = Lagrange{RefHexahedron, order[2]}()^3)
    
        cv = (  left = CellValues(QuadratureRule{RefHexahedron}(4), ip.left, ip.left), # Default geometry interpolation resulting in exception for order=2
                interface = InterfaceCellValues(QuadratureRule{RefQuadrilateral}(4), ip.interface), #, Lagrange{RefQuadrilateral, order}()) # Default geometry interpolation resulting in exception for order=2
                right = CellValues(QuadratureRule{RefHexahedron}(4), ip.right, ip.right))
    
        dh = DofHandler(grid)
        function prepare_sub_dof_handler(key::Symbol)
            s = SubDofHandler(dh, Set{Int}([cellid[key]]))
            if order[1] == order[2] || key == :left
                add!(s, :u, ip[key])
            else
                warntext = "Field :u uses a different interpolation order in another SubDofHandler."
                if key == :interface
                    @test_logs (:warn, warntext) add!(s, :u, ip[key])
                else
                    @test_logs (:warn, warntext) (:warn, warntext) add!(s, :u, ip[key])
                end
            end
            return s
        end
        sdh = ( left = prepare_sub_dof_handler(:left),
                interface = prepare_sub_dof_handler(:interface),
                right = prepare_sub_dof_handler(:right))
        close!(dh)
    
        ch = ConstraintHandler(dh)
        add!(ch, Dirichlet(:u, getfaceset(grid, "∂Ωₗ"), (x,t) -> [0.0, 1.0, 0.0], [1,2,3]))
        add!(ch, Dirichlet(:u, getfaceset(grid, "∂Ωᵣ"), (x,t) -> [2.0, 1.0, 0.0], [1,2,3]))
        close!(ch)
    
        @testset "Accessing shape values etc." begin
            icv = cv[:interface]
            for cc in CellIterator(dh, Set{Int}([cellid[:interface]]))
                reinit!(icv, cc)
            end
            nbf = getnbasefunctions(icv)
            qp = rand(1:getnquadpoints(icv))
            i = rand(1:nbf)
            side, bi = get_side_and_baseindex(icv, i)
            if side == :here
                @test shape_value(icv, qp, i, true) == shape_value(icv.here, qp, bi)
                @test shape_value(icv, qp, i, false) == zero(Vec{3})
                @test shape_value_average(icv, qp, i) == 0.5*shape_value(icv.here, qp, bi)
                @test shape_gradient_average(icv, qp, i) == 0.5*shape_gradient(icv.here, qp, bi)
                @test shape_value_jump(icv, qp, i) == -shape_value(icv.here, qp, bi)
                @test shape_gradient_jump(icv, qp, i) == -shape_gradient(icv.here, qp, bi)
            else
                @test shape_value(icv, qp, i, true) == zero(Vec{3})
                @test shape_value(icv, qp, i, false) == shape_value(icv.there, qp, bi)
                @test shape_value_average(icv, qp, i) == 0.5*shape_value(icv.there, qp, bi)
                @test shape_gradient_average(icv, qp, i) == 0.5*shape_gradient(icv.there, qp, bi)
                @test shape_value_jump(icv, qp, i) == shape_value(icv.there, qp, bi)
                @test shape_gradient_jump(icv, qp, i) == shape_gradient(icv.there, qp, bi)
            end
        end
    
        K = create_sparsity_pattern(dh, ch)
        f = zeros(ndofs(dh))
        assembler = start_assemble(K,f)
    
        for key in (:left, :interface, :right)
            nbf = getnbasefunctions(cv[key])
            Kₑ = zeros(nbf, nbf)
            for cc in CellIterator(dh, Set{Int}([cellid[key]]))
                fill!(Kₑ, 0)
                reinit!(cv[key], cc)
                assemble_test_element!(Kₑ, cv[key])
                assemble!(assembler, celldofs(cc), Kₑ)
            end
        end
    
        apply!(K, f, ch)
        u = K \ f
    
        getind(dofs::Vector) = vcat(collect( collect((i-1)*3+1:(i-1)*3+3) for i in dofs )...)
        tol = 1e-5
    
        for (cellid, leftface, rightface, expected) in ((1,5,3,(0,1,1,0)), (2,1,2,(1,1,1,0)), (3,5,3,(1,2,1,0)))
            cell = getcells(grid, cellid)
            cellfaces = Ferrite.faces(cell)
            cellnodes = Ferrite.get_node_ids(cell)
            dofs = celldofs(dh, cellid)
            leftdofs = dofs[getind(indexin([cellfaces[leftface]...], [cellnodes...]))]
            leftvalues = u[leftdofs]
            rightdofs = dofs[getind(indexin([cellfaces[rightface]...], [cellnodes...]))]
            rightvalues = u[rightdofs]
            @test all(leftvalues[1:3:end-2] .< expected[1]+tol)
            @test all(rightvalues[1:3:end-2] .< expected[2]+tol)
            @test all(vcat(leftvalues[2:3:end-1], rightvalues[2:3:end-1]) .< expected[3]+tol)
            @test all(vcat(leftvalues[3:3:end], rightvalues[3:3:end]) .< expected[4]+tol)
        end
    end
end
    