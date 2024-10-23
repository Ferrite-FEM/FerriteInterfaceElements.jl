using Ferrite, FerriteInterfaceElements, SparseArrays, FerriteGmsh

# grid = togrid("periodic-rve-coarse.msh")
grid = togrid("periodic-rve.msh")

grid = insert_interfaces(grid, ["inclusions", "matrix"]);

ip_bulk = Lagrange{RefTriangle, 1}()
ip_interface = InterfaceCellInterpolation(Lagrange{RefLine, 1}())

qr_bulk = QuadratureRule{RefTriangle}(2)
qr_interface = QuadratureRule{RefLine}(2);

cv_bulk = CellValues(qr_bulk, ip_bulk)
cv_interface = InterfaceCellValues(qr_interface, ip_interface);

dh = DofHandler(grid)
set_bulk = union(getcellset(grid, "inclusions"), getcellset(grid, "matrix"))
set_interface = getcellset(grid, "interfaces")
add!(SubDofHandler(dh, set_bulk),      :u, ip_bulk)
add!(SubDofHandler(dh, set_interface), :u, ip_interface)
close!(dh);

particles = getcellset(grid, "inclusions")
∂Ωᴾ_left   = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "left"))
∂Ωᴾ_right  = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "right"))
∂Ωᴾ_top    = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "top"))
∂Ωᴾ_bottom = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "bottom"));

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ωᴾ_left,   Returns(1.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_right,  Returns(1.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_top,    Returns(0.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_bottom, Returns(0.0)))
close!(ch);

function assemble_element!(Ke::Matrix, cv::CellValues)
    fill!(Ke, 0)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, qp)
        for i in 1:getnbasefunctions(cv)
            ∇δu = shape_gradient(cv, qp, i)
            for j in 1:getnbasefunctions(cv)
                ∇u = shape_gradient(cv, qp, j)
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke
end

function assemble_element!(Ke::Matrix, cv::InterfaceCellValues)
    fill!(Ke, 0)
    for qp in 1:getnquadpoints(cv)
        dΩ = getdetJdV_average(cv, qp)
        for i in 1:getnbasefunctions(cv)
            jump_δu = shape_value_jump(cv, qp, i)
            for j in 1:getnbasefunctions(cv)
                jump_u = shape_value_jump(cv, qp, j)
                Ke[i, j] += (jump_δu * jump_u) * dΩ
            end
        end
    end
    return Ke
end;

function assemble_set!(assembler, set, dh, cv)
    nbf = getnbasefunctions(cv)
    Ke = zeros(nbf, nbf)
    for cc in CellIterator(dh, set)
        reinit!(cv, cc)
        assemble_element!(Ke, cv)
        assemble!(assembler, celldofs(cc), Ke)
    end
    return assembler
end;

function prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
    K = allocate_matrix(dh)
    assembler = start_assemble(K)
    assemble_set!(assembler, set_bulk,      dh, cv_bulk)
    assemble_set!(assembler, set_interface, dh, cv_interface)
    return K, zeros(ndofs(dh))
end;

K, f = prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
apply!(K, f, ch)
u = K \ f;

import Makie, GeometryBasics

function convert_nodes(grid::Grid{dim}) where {dim}
    return collect( GeometryBasics.Point{dim,Float64}(n.x...) for n in grid.nodes )
end

function convert_cells(::Grid)
    bulkcells = filter(c -> !(c isa InterfaceCell), grid.cells) # `InterfaceCell`s are not plotted
    return collect( GeometryBasics.TriangleFace{Int}(cell.nodes...) for cell in bulkcells )
end

function prepare_plotable_mesh(grid::Grid)
    return GeometryBasics.Mesh(convert_nodes(grid), convert_cells(grid))
end

function get_nodal_temperatures(u::AbstractVector, dh::DofHandler)
    grid = dh.grid
    bulkcells = union( getcellset(grid, "inclusions"), getcellset(grid, "matrix") )
    uₙ = zeros(length(grid.nodes))
    for cc in CellIterator(dh, bulkcells)
        cell = getcells(grid, cc.cellid)
        dofs = celldofs(cc)
        nodes = [cell.nodes...]
        uₙ[nodes] .= u[dofs]
    end
    return uₙ
end;

function plot_temperature(u::Vector, dh::DofHandler)
    mesh = prepare_plotable_mesh(dh.grid)
    temperature = get_nodal_temperatures(u, dh)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title="Solution", aspect=Makie.DataAspect())
    p = Makie.mesh!(ax, mesh; color=temperature, colormap=:heat)
    Makie.Colorbar(fig[1,2], p; label="u", labelrotation=0)
    return fig
end;

import CairoMakie
fig = plot_temperature(u, dh)
CairoMakie.save("heat_equation_result.png", fig);

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
