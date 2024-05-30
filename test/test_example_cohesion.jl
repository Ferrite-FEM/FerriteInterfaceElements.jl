function assemble_test_element_cohesion!(Kₑ::Matrix, cv::CellValues, dim::Integer)
    δ(i,j) = i == j ? 1.0 : 0.0
    E = SymmetricTensor{4,dim}( (i,j,k,l) -> 0.5*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k)) ) # Using E=1, ν=0
            
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

function assemble_test_element_cohesion!(Kₑ::Matrix, cv::InterfaceCellValues, ::Integer)
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

@testset "In $(dim)D for order $(order)" for (order, dim, get_grid, get_ip) in (
        (1, 2, prepare_interface_test_grid_2D, prepare_interface_test_vector_interpolation_2D),
        (2, 2, prepare_interface_test_grid_2D, prepare_interface_test_vector_interpolation_2D),
        (1, 3, prepare_interface_test_grid_3D, prepare_interface_test_vector_interpolation_3D),
        (2, 3, prepare_interface_test_grid_3D, prepare_interface_test_vector_interpolation_3D)
        )
    grid = get_grid(order)
    ip, cv = get_ip(order)
    
    dh = DofHandler(grid)
    add!(SubDofHandler(dh, Set{Int}([1])), :u, ip.left)
    add!(SubDofHandler(dh, Set{Int}([2])), :u, ip.interface)
    add!(SubDofHandler(dh, Set{Int}([3])), :u, ip.right)
    close!(dh)
    
    ch = ConstraintHandler(dh)
    if dim == 2
        add!(ch, Dirichlet(:u, getfacetset(grid, "∂Ωₗ"), (x,t) -> [0.0, 1.0], [1,2]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "∂Ωᵣ"), (x,t) -> [2.0, 1.0], [1,2]))
    elseif dim == 3
        add!(ch, Dirichlet(:u, getfacetset(grid, "∂Ωₗ"), (x,t) -> [0.0, 1.0, 0.0], [1,2,3]))
        add!(ch, Dirichlet(:u, getfacetset(grid, "∂Ωᵣ"), (x,t) -> [2.0, 1.0, 0.0], [1,2,3]))
    end
    close!(ch)
    
    K = create_sparsity_pattern(dh, ch)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K,f)
    for (i, key) in enumerate((:left, :interface, :right))
        nbf = getnbasefunctions(cv[key])
        Kₑ = zeros(nbf, nbf)
        for cc in CellIterator(dh, Set{Int}([i]))
            fill!(Kₑ, 0)
            reinit!(cv[key], cc)
            assemble_test_element_cohesion!(Kₑ, cv[key], dim)
            assemble!(assembler, celldofs(cc), Kₑ)
        end
    end
    apply!(K, f, ch)
    u = K \ f

    testinput = dim == 2 ? ((:left,4,2,0,1), (:interface,1,2,1,1), (:right,4,2,1,2)) : ((:left,5,3,0,1), (:interface,1,2,1,1), (:right,5,3,1,2))
    for (i, (key, leftfacet, rightfacet, leftx, rightx)) in enumerate(testinput)
        dofs = celldofs(dh, i)
        facetindices = Ferrite.facetdof_indices(ip[key])
        leftdofs  = dofs[[facetindices[leftfacet]...]]
        rightdofs = dofs[[facetindices[rightfacet]...]]
        
        test_coordinate(given, expected) = abs(given - expected) ≤ 1e-5
            # x coordinates
        @test all( test_coordinate.(u[leftdofs][1:dim:end], [leftx]) )
        @test all( test_coordinate.(u[rightdofs][1:dim:end], [leftx]) )
            # y coordinates
        @test all( test_coordinate.(u[leftdofs][2:dim:end], [1]) )
        @test all( test_coordinate.(u[rightdofs][2:dim:end], [1]) )
        if dim == 3 # z coordinates
            @test all( test_coordinate.(u[leftdofs][3:dim:end], [0]) )
            @test all( test_coordinate.(u[rightdofs][3:dim:end], [0]) )
        end
    end
end