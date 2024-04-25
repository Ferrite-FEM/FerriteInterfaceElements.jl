module FerriteInterfaceElements

using Ferrite
import Ferrite: AbstractCell, AbstractRefShape,
    getnbasefunctions, getngeobasefunctions, getorder, n_components,
    vertexdof_indices, facedof_indices, facedof_interior_indices, celldof_interior_indices,
    vertices, faces, edges,
    adjust_dofs_during_distribution, default_interpolation, default_geometric_interpolation,
    checkbounds, checkquadpoint, function_value_init, function_gradient_init,
    shape_value_type, shape_gradient_type

include("cells.jl")
include("interpolations.jl")
include("cellvalues.jl")

export
    InterfaceCell,
    InterfaceCellInterpolation,
    InterfaceCellValues,
    get_side_and_baseindex,
    shape_value_average,
    shape_gradient_average,
    shape_value_jump,
    shape_gradient_jump,
    function_value_average,
    function_gradient_average,
    function_value_jump,
    function_gradient_jump,
    getdetJdV_average,
    insert_interfaces

end
