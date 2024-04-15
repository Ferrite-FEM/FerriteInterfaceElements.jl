```@meta
EditURL = "../literate-tutorials/heat_equation.jl"
```

# [Heat equation](@id tutorial-heat-equation)

## Introduction

In this example, we solve the stationary heat equation in a domain with particles ($\Omega^\text{P}$)
embedded in matrix ($\Omega^\text{M}$) with a possible temperature jump at the interface ($\Gamma^{^\text{P/M}}$).
For the bulk domains, the following balance law is applied:

```math
 -\nabla \cdot (k \nabla u) = f  \quad \boldsymbol{x} \in \Omega^\text{P} \cup \Omega^\text{M},
```

where $u$ is the unknown temperature field, $k$ the heat conductivity,
$f$ the heat source and $\Omega$ the domain. For simplicity we set $k=1$ and $f = 0$.

For the interface, we formulate a balance law for the normal heat fluxes $q_n^\text{P}$ and $q_n^\text{M}$:

```math
 q_n^\text{P} + q_n^\text{M} = 0  \quad \boldsymbol{x} \in \Gamma^\text{P/M},
```

We define the flux across the interface as:

```math
 q_n := q_n^\text{P} = - k_{if} [\![ u ]\!]
```

where $k_{if}$ is the interface conductivity and $[\![ u ]\!]$ is the jump $u^\text{M} - u^\text{P}$
in temperature from the particle to the matrix domain. For simplicity we set $k_{if} = 1$.

Note that the two sides of an interface are referred to as "here" and "there".
The provided functions use a jump as defined from "here" to "there".
So, in this example the particles are to be considered the side "here".
(Here, this is not that relevant, since the sign of the jump has no influence
on the result in the weak form. In general, this could be relevant.)

The resulting weak form is given given as follows: Find ``u \in \mathbb{U}`` such that
```math
 \int_{\Omega} \nabla \delta u \cdot \nabla u \mathrm{d}\Omega
 \int_{\Gamma^{\text{P}/\text{M}}} [\![ \delta u ]\!] [\![ u ]\!] \mathrm{d}\Gamma
 = 0 \quad \forall \delta u \in \mathbb{T},
```
where $\delta u$ is a test function, and where $\mathbb{U}$ and $\mathbb{T}$ are suitable
trial and test function sets, respectively.

## Commented Program

Now we solve the problem in Ferrite together with FerriteInterfaceElements.
What follows is a program spliced with comments.
The full program, without comments, can be found in the next [section](@ref heat_equation-plain-program).

First we load all the packages we need.

````@example heat_equation
using Ferrite, FerriteInterfaceElements, SparseArrays
Ferrite.facedof_interior_indices(::Lagrange{RefLine}) = Tuple{}() # Hot fix -> TODO: remove ASAP
````

We first load the mesh file [`periodic-rve.msh`](periodic-rve.msh)
([`periodic-rve-coarse.msh`](periodic-rve-coarse.msh) for a coarser mesh). The mesh is
generated with [`gmsh`](https://gmsh.info/), and we read it in as a `Ferrite` grid using
the [`FerriteGmsh`](https://github.com/Ferrite-FEM/FerriteGmsh.jl) package:

````@example heat_equation
using FerriteGmsh
grid = togrid("periodic-rve.msh")
````

````@example heat_equation
grid = redirect_stdout(devnull) do                #hide
    togrid("periodic-rve-coarse.msh") #hide
end                                               #hide
````

So far, this mesh only contains standard cells for the bulk phases.
Thus, we need to insert `InterfaceCell`s at the faces between cells in different phases.
To do so, we can use the function `create_interface_cells!`.
Note that this function does not insert the `InterfaceCell`s into the grid!
(That would require changing the type of `grid.cells`.)
The function adds duplicates of nodes on the interface to the grid, disconnects cells
of different phases and returns cells which can be used to connect them via interfaces.

````@example heat_equation
interface_cells = create_interface_cells!(grid::Grid, "inclusions", "matrix");
nothing #hide
````

Now, we create a new grid with all the needed cells. For convenience, we also add a
new cell set with the interface cells.

````@example heat_equation
n_bulk_cells = length(grid.cells)
n_interface_cells = length(interface_cells)
set_interface = Set{Int}(n_bulk_cells + 1 : n_bulk_cells + n_interface_cells)

grid = Grid(vcat(grid.cells, interface_cells), grid.nodes; cellsets=grid.cellsets, facesets=grid.facesets)
addcellset!(grid, "interface", set_interface);
nothing #hide
````

### Trial and test functions
First, we define an interpolation and a quadrature rule for both bulk and interface cells.

````@example heat_equation
ip_bulk = Lagrange{RefTriangle, 1}()
ip_interface = InterfaceCellInterpolation(Lagrange{RefLine, 1}())

qr_bulk = QuadratureRule{RefTriangle}(2)
qr_interface = QuadratureRule{RefLine}(2);
nothing #hide
````

With these definitions, we can create cell values.

````@example heat_equation
cv_bulk = CellValues(qr_bulk, ip_bulk)
cv_interface = InterfaceCellValues(qr_interface, ip_interface);
nothing #hide
````

### Degrees of freedom
Next we create the `DofHandler`. Here ones needs distinguish the different types of cells.
This can be done by using `SubDofHandlers`.

````@example heat_equation
dh = DofHandler(grid)
set_bulk = union(getcellset(grid, "inclusions"), getcellset(grid, "matrix"))
add!(SubDofHandler(dh, set_bulk),      :u, ip_bulk)
add!(SubDofHandler(dh, set_interface), :u, ip_interface)
close!(dh);
nothing #hide
````

### Boundary conditions
To construct a simple example, but be still be able to observe jumps at the interfaces,
we fix the temperatur to different values on the particle portion of the boundary.

````@example heat_equation
particles = getcellset(grid, "inclusions")
∂Ωᴾ_left   = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "left"))
∂Ωᴾ_right  = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "right"))
∂Ωᴾ_top    = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "top"))
∂Ωᴾ_bottom = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "bottom"));
nothing #hide
````

We set up a `ConstraintHandler` with for fixed values on the four sides.

````@example heat_equation
ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ωᴾ_left,   Returns(0.25)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_right,  Returns(0.75)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_top,    Returns(0.5)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_bottom, Returns(0.0)))
close!(ch);
nothing #hide
````

### Assembling the linear system

Now we can prepare the assembly of the global system of equations.
We start by defining two element assembly functions for the two phases.
To distinguish the phases, we use dispatch und the cell value types.

````@example heat_equation
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
nothing #hide
````

Next, we define a function that performs the assembly for a given set of cells.

````@example heat_equation
function assemble_set!(assembler, set::Set{Int}, dh, cv)
    nbf = getnbasefunctions(cv)
    Ke = zeros(nbf, nbf)
    for cc in CellIterator(dh, set)
        reinit!(cv, cc)
        assemble_element!(Ke, cv)
        assemble!(assembler, celldofs(cc), Ke)
    end
    return assembler
end;
nothing #hide
````

This can then be used in a function for preparing the system of equations.

````@example heat_equation
function prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
    K = create_sparsity_pattern(dh)
    assembler = start_assemble(K)
    assemble_set!(assembler, set_bulk,      dh, cv_bulk)
    assemble_set!(assembler, set_interface, dh, cv_interface)
    return K, zeros(ndofs(dh))
end;
nothing #hide
````

### Solution and visualization
Finally, the system of equations can be assembled and solved.

````@example heat_equation
K, f = prepare_system(dh, cv_bulk, cv_interface, set_bulk, set_interface)
apply!(K, f, ch)
u = K \ f;
nothing #hide
````

To visualize the computed result, we will write a simple function
for plotting the temperature using Makie. The first important step is
translating the Ferrite grid to a representation Makie can use.

````@example heat_equation
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
````

The nex step is to get the nodal temperature values. However, the Ferrite function
`evaluate_at_grid_nodes` does not (yet) work for interface elements.
Thus, we need to rearrange the solution vector ourselves. Althought it is not the
most effiecient way, a simple solution is to iterate over the bulk cells and use the
DOF mapping of each cell.

````@example heat_equation
function get_nodal_temperatures(u::AbstractVector, dh::DofHandler)
    grid = dh.grid
    bulkcells = union( getcellset(grid, "inclusions"), getcellset(grid, "matrix") )
    uₙ = zeros(length(grid.nodes))
    for cc in CellIterator(dh, bulkcells)
        cell = getcells(grid, cc.cellid.x)
        dofs = celldofs(cc)
        nodes = [cell.nodes...]
        uₙ[nodes] .= u[dofs]
    end
    return uₙ
end;
nothing #hide
````

Finally, a function can be defined to create the desired plot.

````@example heat_equation
function plot_temperature(u::Vector, dh::DofHandler)
    mesh = prepare_plotable_mesh(dh.grid)
    temperature = get_nodal_temperatures(u, dh)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1,1]; title="Solution")
    p = Makie.mesh!(ax, mesh; color=temperature, colormap=:heat)
    Makie.Colorbar(fig[1,2], p; label="u", labelrotation=0)
    return fig
end;
nothing #hide
````

Using this function we can create and save a visualization.

````@example heat_equation
import CairoMakie
fig = plot_temperature(u, dh)
CairoMakie.save("heat_equation_result.png", fig);
nothing #hide
````

![](heat_equation_result.png)

## [Plain program](@id heat_equation-plain-program)

Here follows a version of the program without any comments.
The file is also available here: [`heat_equation.jl`](heat_equation.jl).

```julia
using Ferrite, FerriteInterfaceElements, SparseArrays
Ferrite.facedof_interior_indices(::Lagrange{RefLine}) = Tuple{}() # Hot fix -> TODO: remove ASAP

using FerriteGmsh
# grid = togrid("periodic-rve-coarse.msh")
grid = togrid("periodic-rve.msh")

interface_cells = create_interface_cells!(grid::Grid, "inclusions", "matrix");

n_bulk_cells = length(grid.cells)
n_interface_cells = length(interface_cells)
set_interface = Set{Int}(n_bulk_cells + 1 : n_bulk_cells + n_interface_cells)

grid = Grid(vcat(grid.cells, interface_cells), grid.nodes; cellsets=grid.cellsets, facesets=grid.facesets)
addcellset!(grid, "interface", set_interface);

ip_bulk = Lagrange{RefTriangle, 1}()
ip_interface = InterfaceCellInterpolation(Lagrange{RefLine, 1}())

qr_bulk = QuadratureRule{RefTriangle}(2)
qr_interface = QuadratureRule{RefLine}(2);

cv_bulk = CellValues(qr_bulk, ip_bulk)
cv_interface = InterfaceCellValues(qr_interface, ip_interface);

dh = DofHandler(grid)
set_bulk = union(getcellset(grid, "inclusions"), getcellset(grid, "matrix"))
add!(SubDofHandler(dh, set_bulk),      :u, ip_bulk)
add!(SubDofHandler(dh, set_interface), :u, ip_interface)
close!(dh);

particles = getcellset(grid, "inclusions")
∂Ωᴾ_left   = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "left"))
∂Ωᴾ_right  = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "right"))
∂Ωᴾ_top    = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "top"))
∂Ωᴾ_bottom = filter(faceindex -> faceindex[1] in particles, getfaceset(grid, "bottom"));

ch = ConstraintHandler(dh)
add!(ch, Dirichlet(:u, ∂Ωᴾ_left,   Returns(0.25)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_right,  Returns(0.75)))
add!(ch, Dirichlet(:u, ∂Ωᴾ_top,    Returns(0.5)))
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

function assemble_set!(assembler, set::Set{Int}, dh, cv)
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
    K = create_sparsity_pattern(dh)
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
        cell = getcells(grid, cc.cellid.x)
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
    ax = Makie.Axis(fig[1,1]; title="Solution")
    p = Makie.mesh!(ax, mesh; color=temperature, colormap=:heat)
    Makie.Colorbar(fig[1,2], p; label="u", labelrotation=0)
    return fig
end;

import CairoMakie
fig = plot_temperature(u, dh)
CairoMakie.save("heat_equation_result.png", fig);
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

