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
        for func_ip in (fip, fip^3)
            cv = InterfaceCellValues(qr, func_ip, gip)
            @test Ferrite.getngeobasefunctions(cv) == getnbasefunctions(gip)
            @test geometric_interpolation(cv.here) == gip.base
        end
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

@testset "Interface reinitialization allocations" begin
    ip = InterfaceCellInterpolation(Lagrange{RefLine, 2}())
    qr = QuadratureRule{RefLine}(2)
    x = Vec{2,Float64}.([(-1, 0), (1, 0), (-1, 0), (1, 0), (0, 0), (0, 0)])
    for shared in (true, false)
        cv = InterfaceCellValues(qr, ip; use_same_cv=shared)
        reinit!(cv, x)
        @test (@allocated reinit!(cv, x)) == 0
    end
end

@testset "Mixed solution and geometry orders" begin
    for (shape, cells, dim) in ((RefLine, (Line, QuadraticLine), 2),
                               (RefTriangle, (Triangle, QuadraticTriangle), 3),
                               (RefQuadrilateral, (Quadrilateral, QuadraticQuadrilateral), 3))
        for forder in (1, 2), gorder in (1, 2), vectorized in (false, true), shared in (false, true)
            base_fip = Lagrange{shape,forder}()
            base_gip = Lagrange{shape,gorder}()
            fip = InterfaceCellInterpolation(base_fip)
            gip = InterfaceCellInterpolation(base_gip)
            qr = QuadratureRule{shape}(2)
            cv = InterfaceCellValues(qr, vectorized ? fip^dim : fip, gip;
                                     use_same_cv=shared, include_R=true)
            xh = [Vec{dim}(i -> i < dim ? ξ[i] : 0.0) for ξ in Ferrite.reference_coordinates(base_gip)]
            xt = [2x + Vec{dim}(i -> i == dim ? 1.0 : 0.0) for x in xh]
            n = length(xh)
            C = cells[gorder]
            cell = InterfaceCell(C(Tuple(1:n)), C(Tuple(n+1:2n)))
            x = vcat(xh, xt)[collect(cell.nodes)]
            reinit!(cv, x)
            for (actual, coords) in ((cv.here, xh), (cv.there, shared ? xh : xt))
                expected = CellValues(qr, vectorized ? base_fip^dim : base_fip, base_gip^dim)
                reinit!(expected, coords)
                for qp in 1:getnquadpoints(cv)
                    @test getdetJdV(actual, qp) ≈ getdetJdV(expected, qp)
                    @test midplane_rotation(cv, qp) ≈ one(Tensor{2,dim,Float64})
                    for i in 1:getnbasefunctions(actual)
                        @test shape_gradient(actual, qp, i) ≈ shape_gradient(expected, qp, i)
                    end
                end
            end
            @test_throws ArgumentError reinit!(cv, x[1:end-1])
            @test_throws ArgumentError reinit!(cv, vcat(x, x[1:1]))
        end
    end
end

@testset "Invalid interface solution lengths" begin
    ip = InterfaceCellInterpolation(Lagrange{RefLine, 1}())
    cv = InterfaceCellValues(QuadratureRule{RefLine}(2), ip)
    for f in (function_value, function_gradient), here in (true, false)
        @test_throws ArgumentError f(cv, 1, ones(3), here)
        @test_throws ArgumentError f(cv, 1, ones(5), here)
        @test_throws ArgumentError f(cv, 1, ones(8), here, 2:4)
    end
end

@testset "Independent side cell values" begin
    ip = InterfaceCellInterpolation(Lagrange{RefLine, 1}())
    qr = QuadratureRule{RefLine}(2)
    for fip in (ip, ip^2)
        cv = InterfaceCellValues(qr, fip; use_same_cv=false)
        x = Vec{2,Float64}.([(-1, 0), (1, 0), (-2, 1), (2, 1)])
        reinit!(cv, x)
        @test cv.here !== cv.there
        for qp in 1:getnquadpoints(cv)
            @test getdetJdV(cv.there, qp) ≈ 2getdetJdV(cv.here, qp)
            for i in 1:getnbasefunctions(cv.here)
                @test shape_gradient(cv.here, qp, i) ≈ 2shape_gradient(cv.there, qp, i)
            end
        end
        reinit!(cv.there, 3x[3:4])
        @test getdetJdV(cv.here, 1) ≈ 1.0
        @test getdetJdV(cv.there, 1) ≈ 6.0
    end
end

function test_side_evaluation(cv, x)
    reinit!(cv, x)
    n = getnbasefunctions(cv)
    u = sin.(collect(1:n))
    for (f, shape) in ((function_value, shape_value), (function_gradient, shape_gradient)), here in (true, false)
        for qp in 1:getnquadpoints(cv)
            expected = sum(shape(cv, qp, i, here) * u[i] for i in 1:n)
            @test f(cv, qp, u, here) ≈ expected
            @test f(cv, qp, vcat(99.0, u, 99.0), here, 2:n+1) ≈ expected
            @test f(cv, qp, reverse(u), here, collect(n:-1:1)) ≈ expected
        end
        @test_throws ErrorException f(cv, 0, u, here)
        @test_throws ErrorException f(cv, getnquadpoints(cv)+1, u, here)
        @test_throws BoundsError f(cv, 1, u, here, 2:n+1)
        f(cv, 1, u, here)
        @test (@allocated f(cv, 1, u, here)) == 0
    end
end

@testset "Delegated interface field evaluation" begin
    for (shape, dim) in ((RefLine, 2), (RefTriangle, 3), (RefQuadrilateral, 3))
        gip = InterfaceCellInterpolation(Lagrange{shape,1}())
        xh = [Vec{dim}(i -> i < dim ? ξ[i] : 0.0) for ξ in Ferrite.reference_coordinates(gip.base)]
        xt = [2x + Vec{dim}(i -> i == dim ? 1.0 : 0.0) for x in xh]
        for order in (1, 2), vdim in (1, dim), shared in (false, true)
            ip = InterfaceCellInterpolation(Lagrange{shape,order}())
            cv = InterfaceCellValues(QuadratureRule{shape}(2), vdim == 1 ? ip : ip^vdim, gip;
                                     use_same_cv=shared)
            test_side_evaluation(cv, vcat(xh, xt))
        end
    end
end
