function assemble_test_element_diffusion!(Kₑ::Matrix, cv::CellValues)
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

function assemble_test_element_diffusion!(Kₑ::Matrix, cv::InterfaceCellValues)
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

@testset "In $(dim)D for order $(order)" for (order, dim, get_grid, get_ip) in (
        (1, 2, prepare_interface_test_grid_2D, prepare_interface_test_scalar_interpolation_2D),
        (2, 2, prepare_interface_test_grid_2D, prepare_interface_test_scalar_interpolation_2D),
        (1, 3, prepare_interface_test_grid_3D, prepare_interface_test_scalar_interpolation_3D),
        (2, 3, prepare_interface_test_grid_3D, prepare_interface_test_scalar_interpolation_3D)
        )
    grid = get_grid(order)
    ip, cv = get_ip(order)

    dh = DofHandler(grid)
    add!(SubDofHandler(dh, Set{Int}([1])), :c, ip.left)
    add!(SubDofHandler(dh, Set{Int}([2])), :c, ip.interface)
    add!(SubDofHandler(dh, Set{Int}([3])), :c, ip.right)
    close!(dh)
    
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:c, getfacetset(grid, "∂Ωₗ"), (x,t) -> 3.0))
    add!(ch, Dirichlet(:c, getfacetset(grid, "∂Ωᵣ"), (x,t) -> 0.0))
    close!(ch)
    
    K = allocate_matrix(dh, ch)
    f = zeros(ndofs(dh))
    assembler = start_assemble(K,f)
    for (i, key) in enumerate((:left, :interface, :right))
        nbf = getnbasefunctions(cv[key])
        Kₑ = zeros(nbf, nbf)
        for cc in CellIterator(dh, Set{Int}([i]))
            fill!(Kₑ, 0)
            reinit!(cv[key], cc)
            assemble_test_element_diffusion!(Kₑ, cv[key])
            assemble!(assembler, celldofs(cc), Kₑ)
        end
    end
    apply!(K, f, ch)
    u = K \ f
    
    expectedvalue = 3.0
    testinput = dim == 2 ? ((:left,4,2), (:interface,1,2), (:right,4,2)) : ((:left,5,3), (:interface,1,2), (:right,5,3))
    for (i, (key, leftfacet, rightfacet)) in enumerate(testinput)
        dofs = celldofs(dh, i)
        facetindices = Ferrite.facetdof_indices(ip[key])
        leftdofs  = dofs[[facetindices[leftfacet]...]]
        rightdofs = dofs[[facetindices[rightfacet]...]]
        
        @test all(u[leftdofs]  .≈ expectedvalue)
        @test all(u[rightdofs] .≈ expectedvalue-1)
        expectedvalue -= 1
    end
end