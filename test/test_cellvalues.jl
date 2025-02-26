@testset "InterfaceCellValues" begin
    qr = 
    for (refshape, dim, iporder, giporder) in ( 
            (RefTriangle, 3, 1, 1), (RefTriangle, 3, 1, 2), (RefTriangle, 3, 2, 1), (RefTriangle, 3, 2, 2),
            (RefLine, 2, 1, 1), (RefLine, 2, 1, 2), (RefLine, 2, 2, 1), (RefLine, 2, 2, 2))
        qr  = QuadratureRule{refshape}(2)
        fip = InterfaceCellInterpolation(Lagrange{refshape, iporder}())
        gip = InterfaceCellInterpolation(Lagrange{refshape, giporder}())
        @test InterfaceCellValues(qr, fip, gip^dim) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip^dim, gip^dim) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip^dim, gip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip, gip) isa InterfaceCellValues
        @test InterfaceCellValues(Float32, qr, fip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip) isa InterfaceCellValues
        @test InterfaceCellValues(qr, fip; use_same_cv=false) isa InterfaceCellValues
    end

    for (refshape, dim) in ((RefTriangle,3), (RefLine,2))
        qr  = QuadratureRule{refshape}(2)
        ip = InterfaceCellInterpolation(Lagrange{refshape, 1}())
        cv = InterfaceCellValues(qr, ip)

        @test getnbasefunctions(cv) == 2*dim
        @test Ferrite.getngeobasefunctions(cv) == 2*dim

        x = repeat([rand(Vec{dim}) for _ in 1:dim], 2)
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
    ip = InterfaceCellInterpolation(Lagrange{RefLine, 1}())
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