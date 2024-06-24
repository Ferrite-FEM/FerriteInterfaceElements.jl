@testset "InterfaceCellValues" begin
    qr = QuadratureRule{RefTriangle}(1)
    for (fip, gip) in ( (InterfaceCellInterpolation(Lagrange{RefTriangle, 1}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 2}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 1}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 2}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 2}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 2}())),
                      )
        @test InterfaceCellValues(qr, fip, gip^3) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip^3, gip^3) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip^3, gip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip, gip) isa InterfaceCellValues
        @test InterfaceCellValues(Float32, qr, fip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip; use_same_cv=false) isa InterfaceCellValues
    end

    ip = InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())
    cv = InterfaceCellValues(qr, ip)

    @test getnbasefunctions(cv) == 6
    @test Ferrite.getngeobasefunctions(cv) == 6

    x = repeat([rand(Vec{3}), rand(Vec{3}), rand(Vec{3})], 2)
    reinit!(cv, x)
    nbf = getnbasefunctions(cv)
    here, there = rand(2)
    u = vcat(ones(nbf÷2).*here, ones(nbf÷2).*there)
    for qp in 1:getnquadpoints(cv)
        @test function_value(cv, qp, u; here=true) ≈ here
        @test function_value(cv, qp, u; here=false) ≈ there
        @test all(abs.(function_gradient(cv, qp, u; here=true)) .≤ 1e-14)
        @test all(abs.(function_gradient(cv, qp, u; here=false)) .≤ 1e-14)
        @test function_value_average(cv, qp, u) ≈ (here + there)/2
        @test function_value_jump(cv, qp, u) ≈ here - there
        @test all(abs.(function_gradient_average(cv, qp, u)) .≤ 1e-14)
        @test all(abs.(function_gradient_jump(cv, qp, u)) .≤ 1e-14)
        @test getdetJdV_average(cv, qp) == (getdetJdV(cv.here, qp) + getdetJdV(cv.there, qp)) / 2
    end
end

@testset "Consistency with InterfaceValues" begin
    qr = QuadratureRule{RefTriangle}(1)
    ip = InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())
    cv = InterfaceCellValues(qr, ip)

    ip = DiscontinuousLagrange{RefTetrahedron, 1}()
    qr = FacetQuadratureRule{RefTetrahedron}(1)
    iv = InterfaceValues(qr, ip)

    x₁ = [ Vec((0.0, 0.0, 0.0)), Vec((rand(), 0.0, 0.0)), Vec((0.0, rand(), 0.0)), Vec((0.0, 0.0, rand())) ]
    x₂ = [ Vec((1.0, 1.0, 1.0)), x₁[2], x₁[3], x₁[4] ]
    xᵢ = repeat([x₁[2], x₁[3], x₁[4]], 2)
    u₁, u₂ = rand(4), rand(4)
    uᵢ = [ u₁[2], u₁[3], u₁[4], u₂[2], u₂[3], u₂[4]]

    c₁ = Tetrahedron((1,2,3,4))
    c₂ = Tetrahedron((5,2,3,4))

    reinit!(cv, xᵢ)
    reinit!(iv, c₁, x₁, 3, c₂, x₂, 3)

    for qp in 1:getnquadpoints(cv)
        #@test getdetJdV(cv, qp) == getdetJdV(iv, qp) <- TODO
        @test function_value(cv, qp, uᵢ; here=true)  == function_value(iv, qp, vcat(u₁, u₂); here=true)
        @test function_value(cv, qp, uᵢ; here=false) == function_value(iv, qp, vcat(u₁, u₂); here=false)
        @test function_value_jump(cv, qp, uᵢ) == function_value_jump(iv, qp, vcat(u₁, u₂))
        #@test function_gradient_jump(cv, qp, uᵢ) == function_gradient_jump(iv, qp, vcat(u₁, u₂)) <- TODO
        #@test getnormal(cv, qp) == getnormal(iv, qp) <- TODO
    end
end