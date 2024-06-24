"""
    InterfaceCell(here::AbstractCell, there::AbstractCell) <: AbstractCell

An `InterfaceCell` is a cell based on two cells of lower dimension representing the two facets.
The two base cells need to be of the same type and the order of nodes needs to match, e.g.:
```
1---2 "here"
4---3 "there"
InterfaceCell(Line((1,2)), Line((4,3)))
```

# Fields
- `here::AbstractCell`: cell representing the facet "here"
- `there::AbstractCell`: cell representing the facet "there"
- `nodes`::NTuple: tuple with all node indices in appropriate order: vertex nodes "here", vertex nodes "there", facet nodes "here", ...
"""
struct InterfaceCell{shape, C, N} <: AbstractCell{shape}
    here::C
    there::C
    nodes::NTuple{N,Int}

    function InterfaceCell{shape, C}(here::C, there::C) where {shape<:AbstractRefShape, C<:AbstractCell}
        sni = get_sides_and_base_indices(C)
        nodes = ntuple( i -> sni[i][1] == :here ? here.nodes[sni[i][2]] : there.nodes[sni[i][2]], length(sni))
        new{shape, C, length(nodes)}(here, there, nodes)
    end
end

function InterfaceCell(here::C, there::C) where {baseshape<:AbstractRefShape, C<:AbstractCell{baseshape}}
    shape = get_interface_cell_shape(baseshape)
    return InterfaceCell{shape, C}(here, there)
end

"""
    get_interface_cell_shape(::Type{<:AbstractRefShape})

Return the shape of an interface given a base reference shape.
E.g. given `RefTriangle`, `RefPrism` is returned, meaning two triangles form an interface based on a prism.
"""
get_interface_cell_shape(::Type{RefLine}) = RefQuadrilateral
get_interface_cell_shape(::Type{RefTriangle}) = RefPrism
get_interface_cell_shape(::Type{RefQuadrilateral}) = RefHexahedron

"""
    get_sides_and_base_indices(::InterfaceCell)
    get_sides_and_base_indices(::AbstractCell)
    get_sides_and_base_indices(::Type{<:AbstractCell})

Return a tuple containing tuples of a symbol (:here or :there) and an integer.
The index of the outer tuple represents the node index.
In the inner tuple, the symbol represents the side the node is on 
and the integer represents the nodes index in the base cell.
"""
get_sides_and_base_indices(::InterfaceCell{<:Any, C}) where {C<: AbstractCell} = get_sides_and_base_indices(C)
get_sides_and_base_indices(::C) where {C <: AbstractCell} = get_sides_and_base_indices(C)

get_sides_and_base_indices(::Type{Line}) = ((:here,1), (:here,2), (:there,1), (:there,2))
get_sides_and_base_indices(::Type{QuadraticLine}) = ((:here,1), (:here,2), (:there,1), (:there,2), (:here,3), (:there,3))

get_sides_and_base_indices(::Type{Triangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3))
get_sides_and_base_indices(::Type{QuadraticTriangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3), (:here,4), (:here,5), (:here,6), (:there,4), (:there,5), (:there,6))

get_sides_and_base_indices(::Type{Quadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4))
get_sides_and_base_indices(::Type{QuadraticQuadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4), (:here,5), (:here,6), (:here,7), (:here,8), (:there,5), (:there,6), (:there,7), (:there,8), (:here,9), (:there,9))

Ferrite.vertices(c::InterfaceCell) = (vertices(c.here)..., vertices(c.there)...)
Ferrite.edges(c::InterfaceCell) = (edges(c.here)..., edges(c.there)...)
Ferrite.faces(c::InterfaceCell{shape}) where {shape<:Ferrite.AbstractRefShape{3}} = (vertices(c.here), vertices(c.there))
Ferrite.facets(c::InterfaceCell) = (vertices(c.here), vertices(c.there))
