@testset "InterfaceCellValues" begin
    qr = QuadratureRule{RefTriangle}(1)
    for (fip, gip) in ( (InterfaceCellInterpolation(Lagrange{RefTriangle, 1}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 2}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 1}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 2}())),
                        (InterfaceCellInterpolation(Lagrange{RefTriangle, 2}()), InterfaceCellInterpolation(Lagrange{RefTriangle, 2}())),
                      )
        @inferred InterfaceCellValues(qr, fip, gip^3)
        #@inferred InterfaceCellValues(qr, fip^3, gip^3)
        #@inferred InterfaceCellValues(qr, fip^3, gip)
        @inferred InterfaceCellValues(qr, fip, gip)
        @inferred InterfaceCellValues(Float32, qr, fip)
        @inferred InterfaceCellValues(qr, fip)
        @inferred InterfaceCellValues(qr, fip; use_same_cv=Val(false))
        @inferred InterfaceCellValues(qr, fip; use_same_cv=Val(true))
        @inferred InterfaceCellValues(qr, fip; include_R=Val(false))
        @inferred InterfaceCellValues(qr, fip; include_R=Val(true))
    end

    ip = InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())
    for cv in (InterfaceCellValues(qr, ip), InterfaceCellValues(qr, ip; use_same_cv=Val(false)), InterfaceCellValues(qr, ip; include_R=Val(true)))
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
            n = @allocated function_value_jump(cv, qp, u)
            @test n == 0
        end
    end

    qr  = QuadratureRule{RefTriangle}(2)
    ip = InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())
    cv = InterfaceCellValues(qr, ip; include_R=true)

    x  = repeat([Vec{3}((1.0,0.0,0.0)), Vec{3}((0.0,1.0,0.0)), Vec{3}((0.0,0.0,1.0))], 2)
    n  = Vec{3}(( sqrt(1/3),  sqrt(1/3),   sqrt(1/3)))
    t₁ = Vec{3}(( sqrt(1/2),  0.0,        -sqrt(1/2)))
    t₂ = Vec{3}((-sqrt(1/6), 2*sqrt(1/6), -sqrt(1/6)))
    reinit!(cv, x)
    for qp in 1:getnquadpoints(cv)
        R = midplane_rotation(cv, qp)
        @test tdot(R) ≈ one(R) # Fails e.g. when dx/dξ₁ not perpendicular to dx/ξ₂
        @test R⋅Vec{3}((1.0,0.0,0.0)) ≈ t₁
        @test R⋅Vec{3}((0.0,1.0,0.0)) ≈ t₂
        @test R⋅Vec{3}((0.0,0.0,1.0)) ≈ n
    end

    qr  = QuadratureRule{RefLine}(2)
    for ip in (InterfaceCellInterpolation(Lagrange{RefLine, 1}()), InterfaceCellInterpolation(Lagrange{RefLine, 1}())^2)
        cv = InterfaceCellValues(qr, ip; include_R=true)
        x  = repeat([Vec{2}((1.0,0.0)), Vec{2}((0.0,1.0))], 2)
        n  = Vec{2}((-sqrt(1/2), -sqrt(1/2)))
        t  = Vec{2}((-sqrt(1/2),  sqrt(1/2)))
        reinit!(cv, x)
        for qp in 1:getnquadpoints(cv)
            R = midplane_rotation(cv, qp)
            @test tdot(R) ≈ one(R) # Fails e.g. when dx/dξ₁ not perpendicular to dx/ξ₂
            @test R⋅Vec{2}((1.0,0.0)) ≈ t
            @test R⋅Vec{2}((0.0,1.0)) ≈ n
        end
    end
end

@testset "reinit! with cell argument" begin
    qr = QuadratureRule{RefTriangle}(1)
    ip = InterfaceCellInterpolation(Lagrange{RefTriangle, 1}())
    cell = InterfaceCell(Triangle((1, 2, 3)), Triangle((4, 5, 6)))
    x = repeat([rand(Vec{3}), rand(Vec{3}), rand(Vec{3})], 2)

    for kwargs in ((), (; include_R = Val(true)), (; use_same_cv = Val(false)))
        cv_2arg = InterfaceCellValues(qr, ip; kwargs...)
        cv_3arg = InterfaceCellValues(qr, ip; kwargs...)
        cv_nothing = InterfaceCellValues(qr, ip; kwargs...)

        reinit!(cv_2arg, x)
        reinit!(cv_3arg, cell, x)
        reinit!(cv_nothing, nothing, x)

        for qp in 1:getnquadpoints(cv_2arg)
            @test getdetJdV_average(cv_3arg, qp)   == getdetJdV_average(cv_2arg, qp)
            @test getdetJdV_average(cv_nothing, qp) == getdetJdV_average(cv_2arg, qp)
            for i in 1:getnbasefunctions(cv_2arg)
                @test shape_value_jump(cv_3arg, qp, i)    == shape_value_jump(cv_2arg, qp, i)
                @test shape_value_jump(cv_nothing, qp, i) == shape_value_jump(cv_2arg, qp, i)
            end
        end
    end
end

@testset "reinit! with superparametric field (ip order 2, geometry order 1)" begin
    # Quadratic field, linear geometry: base_indices_* (function space, 3 per side)
    # and the coordinate vector (4 geometric nodes) have different lengths.
    ip     = InterfaceCellInterpolation(Lagrange{RefLine, 2}())
    ip_geo = VectorizedInterpolation{2}(InterfaceCellInterpolation(Lagrange{RefLine, 1}()))
    qr     = QuadratureRule{RefLine}(2)

    p1 = Vec{2}((0.0, 0.0))
    p2 = Vec{2}((2.0, 1.0))
    x  = [p1, p2, p1, p2]          # node order: (here,1), (here,2), (there,1), (there,2)

    t = (p2 - p1) / norm(p2 - p1)
    n = Vec{2}((-t[2], t[1]))

    @testset "include_R = false" begin
        cv = InterfaceCellValues(Float64, qr, ip, ip_geo; include_R = Val(false))
        @test length(cv.base_indices_here) == 3
        @test length(x) == 4
        reinit!(cv, x)
        for qp in 1:getnquadpoints(cv.here)
            @test getdetJdV_average(cv, qp) > 0
        end
    end

    @testset "include_R = true" begin
        cv = InterfaceCellValues(Float64, qr, ip, ip_geo; include_R = Val(true))
        reinit!(cv, x)
        for qp in 1:getnquadpoints(cv.here)
            R = midplane_rotation(cv, qp)
            @test R⋅Vec{2}((1.0,0.0)) ≈ t
            @test R⋅Vec{2}((0.0,1.0)) ≈ n
            @test det(R) ≈ 1.0
            @test R' ⋅ R ≈ one(Tensor{2, 2})
        end
    end
end