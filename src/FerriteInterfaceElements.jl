module FerriteInterfaceElements

import Ferrite: Ferrite, AbstractCell, AbstractRefShape, AbstractCellValues,
    RefLine, RefQuadrilateral, RefTriangle, RefPrism, RefHexahedron,
    Line, QuadraticLine, Triangle, QuadraticTriangle, Quadrilateral, QuadraticQuadrilateral, Tetrahedron, Hexahedron,
    ScalarInterpolation, VectorizedInterpolation, Lagrange,
    CellValues, QuadratureRule, CellCache, Grid, ExclusiveTopology, FacetIndex, Vec,
    getnbasefunctions, getngeobasefunctions, getorder, n_components, getrefshape,
    vertexdof_indices, edgedof_interior_indices, facedof_interior_indices, volumedof_interior_indices,
    vertices, facets, edges,
    adjust_dofs_during_distribution, default_geometric_interpolation, geometric_interpolation,
    checkbounds, checkquadpoint, function_value_init, function_gradient_init,
    shape_value_type, shape_gradient_type,
    reinit!, getnquadpoints, getdetJdV, shape_value, shape_gradient, function_value, function_gradient,
    getcellset, getcells, getneighborhood,
    ×, norm
import OrderedCollections: OrderedSet
import StaticArraysCore: SMatrix, SVector

include("cells.jl")
export InterfaceCell

include("interpolations.jl")
export InterfaceCellInterpolation

include("cellvalues.jl")
export InterfaceCellValues,
    get_side_and_baseindex,
    shape_value_average, shape_gradient_average,
    shape_value_jump, shape_gradient_jump,
    function_value_average, function_gradient_average,
    function_value_jump, function_gradient_jump,
    midplane_rotation,
    getdetJdV_average

include("grid.jl")
export insert_interfaces

end
