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
        @test function_value(cv, qp, u, true) ≈ here
        @test function_value(cv, qp, u, false) ≈ there
        @test all(abs.(function_gradient(cv, qp, u, true)) .≤ 1e-14)
        @test all(abs.(function_gradient(cv, qp, u, false)) .≤ 1e-14)
        @test function_value_average(cv, qp, u) ≈ (here + there)/2
        @test function_value_jump(cv, qp, u) ≈ there - here
        @test all(abs.(function_gradient_average(cv, qp, u)) .≤ 1e-14)
        @test all(abs.(function_gradient_jump(cv, qp, u)) .≤ 1e-14)
        @test getdetJdV_average(cv, qp) == (getdetJdV(cv.here, qp) + getdetJdV(cv.there, qp)) / 2
    end

    cv = InterfaceCellValues(qr, ip; include_R=true)
    x = repeat([Vec{3}((0.0,0.0,0.0)), Vec{3}((1.0,0.0,0.0)), Vec{3}((0.0,1.0,0.0))], 2)
    reinit!(cv, x)
    for qp in 1:getnquadpoints(cv)
        @test midplane_rotation(cv, qp) ≈ [ 0.0  sqrt(2)/2 0.0;
                                           -1.0 -sqrt(2)/2 0.0;
                                            0.0  0.0       1.0]
    end
end