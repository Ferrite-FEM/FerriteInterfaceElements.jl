# # [Heat equation](@id tutorial-heat-equation)
#
# ## Introduction
#
# In this example, we solve the stationary heat equation in a domain with particles ($\Omega^\text{P}$)
# embedded in matrix ($\Omega^\text{M}$) with a possible temperature jump at the interface ($\Gamma^{^\text{P/M}}$).
# For the bulk domains, the following balance law is applied:
#
# ```math
#  -\nabla \cdot (k \nabla u) = f  \quad \boldsymbol{x} \in \Omega^\text{P} \cup \Omega^\text{M},
# ```
#
# where $u$ is the unknown temperature field, $k$ the heat conductivity,
# $f$ the heat source and $\Omega$ the domain. For simplicity we set $k=1$ and $f = 0$.
# 
# For the interface, we formulate a balance law for the normal heat fluxes $q_n^\text{P}$ and $q_n^\text{M}$:
#
# ```math
#  q_n^\text{P} + q_n^\text{M} = 0  \quad \boldsymbol{x} \in \Gamma^\text{P/M},
# ```
#
# We define the flux across the interface as:
#
# ```math
#  q_n := q_n^\text{P} = - k_{if} [\![ u ]\!]
# ```
#
# where $k_{if}$ is the interface conductivity and $[\![ u ]\!]$ is the jump $u^\text{M} - u^\text{P}$
# in temperature from the particle to the matrix domain. For simplicity we set $k_{if} = 1$.
#
# Note that the two sides of an interface are referred to as "here" and "there".
# The provided functions use a jump as defined from "here" to "there".
# So, in this example the particles are to be considered the side "here".
# (Here, this is not really relevant, since the sign of the jump has no influence
# on the result in the weak form. In general, this could be relevant.)
#
# The resulting weak form is given given as follows: Find ``u \in \mathbb{U}`` such that
# ```math
#  \int_{\Omega} \nabla \delta u \cdot \nabla u \mathrm{d}\Omega
#  + \int_{\Gamma^{\text{P}/\text{M}}} [\![ \delta u ]\!] [\![ u ]\!] \mathrm{d}\Gamma
#  = 0 \quad \forall \delta u \in \mathbb{T},
# ```
# where $\delta u$ is a test function, and where $\mathbb{U}$ and $\mathbb{T}$ are suitable
# trial and test function spaces, respectively.
#-
# ## Commented Program
#
# Now we solve the problem in Ferrite together with FerriteInterfaceElements.
# What follows is a program spliced with comments.
#md # The full program, without comments, can be found in the next [section](@ref heat_equation-plain-program).
#
# First we load all the packages we need for the computation.

using Ferrite, FerriteInterfaceElements, SparseArrays, FerriteGmsh

# Then, we load the mesh file [`periodic-rve.msh`](periodic-rve.msh)
# ([`periodic-rve-coarse.msh`](periodic-rve-coarse.msh) for a coarser mesh). The mesh is
# generated with [`gmsh`](https://gmsh.info/), and we read it in as a `Ferrite` grid using
# the [`FerriteGmsh`](https://github.com/Ferrite-FEM/FerriteGmsh.jl) package:

#src notebook: use coarse mesh to decrease build time
#src   script: use the fine mesh
#src markdown: use the coarse mesh to decrease build time, but make it look like the fine
#nb ## grid = togrid("periodic-rve.msh")
#nb grid = togrid("periodic-rve-coarse.msh")
#jl ## grid = togrid("periodic-rve-coarse.msh")
#jl grid = togrid("periodic-rve.msh")
#md grid = togrid("periodic-rve.msh")
#-
#md grid = redirect_stdout(devnull) do                #hide
#md     togrid("periodic-rve-coarse.msh") #hide
#md end                                               #hide

grid = togrid("periodic-rve.msh") #src

# So far, this mesh only contains standard cells for the bulk phases.
# Thus, we need to insert `InterfaceCell`s at the faces between cells in different phases.
# To do so, we can use the function `insert_interfaces`.
# This function creates a new grid with interfaces between subdomains which are defined by
# names of cellsets which are passed as arguments.
# The resulting grid includes new cellsets: `"interfaces"` for all interfaces and
# `domain1-domain2-interface"` for interfaces between each pair of subdomains.
# Note that original cell and face sets are preserved, however, node sets are not.

grid = insert_interfaces(grid, ["inclusions", "matrix"]);

# ### Trial and test functions
# First, we define an interpolation and a quadrature rule for both bulk and interface cells.

ip_bulk = Lagrange{RefTriangle, 1}()
ip_interface = InterfaceCellInterpolation(Lagrange{RefLine, 1}())

qr_bulk = QuadratureRule{RefTriangle}(2)
qr_interface = QuadratureRule{RefLine}(2);

# With these definitions, we can create cell values.

cv_bulk = CellValues(qr_bulk, ip_bulk)
cv_interface = InterfaceCellValues(qr_interface, ip_interface);

# ### Degrees of freedom
# Next we create the `DofHandler`. Here, one needs distinguish the different types of cells.
# This can be done by using `SubDofHandlers`.

dh = DofHandler(grid)
set_bulk = union(getcellset(grid, "inclusions"), getcellset(grid, "matrix"))
set_interface = getcellset(grid, "interfaces")
add!(SubDofHandler(dh, set_bulk),      :u, ip_bulk)
add!(SubDofHandler(dh, set_interface), :u, ip_interface)
close!(dh);

# ### Boundary conditions
# To construct a simple example and still be able to observe jumps at the interfaces,
# we fix the temperatur to different values on the particle portion of the boundary.

particles = getcellset(grid, "inclusions")
∂Ωᴾ_left   = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "left"))
∂Ωᴾ_right  = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "right"))
∂Ωᴾ_top    = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "top"))
∂Ωᴾ_bottom = filter(facetindex -> facetindex[1] in particles, getfacetset(grid, "bottom"));

# We set up a `ConstraintHandler` with fixed values on the four sides.

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ωᴾ_left,   Returns(1.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_right,  Returns(1.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_top,    Returns(0.0)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_bottom, Returns(0.0)))
close!(ch);

# ### Assembling the linear system
#
# Now we can prepare the assembly of the global system of equations.
# We start by defining two element assembly functions for the two phases.
# To distinguish the phases, we use dispatch on the cell value types.

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

# Next, we define a function that performs the assembly for a given set of cells.

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

# This can then be used in a function for preparing the system of equations.

function prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
    K = allocate_matrix(dh)
    assembler = start_assemble(K)
    assemble_set!(assembler, set_bulk,      dh, cv_bulk)
    assemble_set!(assembler, set_interface, dh, cv_interface)
    return K, zeros(ndofs(dh))
end;

# ### Solution and visualization
# Finally, the system of equations can be assembled and solved.

K, f = prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
apply!(K, f, ch)
u = K \ f;

# To visualize the computed result, we will write a simple function
# for plotting the temperature using Makie. The first important step is
# translating the Ferrite grid to a representation Makie can use.
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

# The nex step is to get the nodal temperature values. However, the Ferrite function
# `evaluate_at_grid_nodes` does not (yet) work for interface elements.
# Thus, we need to rearrange the solution vector ourselves. Althought it is not the
# most efficient way, a simple solution is to iterate over the bulk cells and use the 
# DOF mapping of each cell.

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

# Finally, a function can be defined to create the desired plot.

function plot_temperature(u::Vector, dh::DofHandler)
    mesh = prepare_plotable_mesh(dh.grid)
    temperature = get_nodal_temperatures(u, dh)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title="Solution", aspect=Makie.DataAspect())
    p = Makie.mesh!(ax, mesh; color=temperature, colormap=:heat)
    Makie.Colorbar(fig[1,2], p; label="u", labelrotation=0)
    return fig
end;

# Using this function we can create and save a visualization.

import CairoMakie
fig = plot_temperature(u, dh)
CairoMakie.save("heat_equation_result.png", fig);

# ![](heat_equation_result.png)

#md # ## [Plain program](@id heat_equation-plain-program)
#md #
#md # Here follows a version of the program without any comments.
#md # The file is also available here: [`heat_equation.jl`](heat_equation.jl).
#md #
#md # ```julia
#md # @__CODE__
#md # ```
