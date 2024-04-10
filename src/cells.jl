"""
    InterfaceCell(here::AbstractCell, there::AbstractCell) <: AbstractCell

An `InterfaceCell` is a cell based on two cells of lower dimension representing the two faces.
The two base cells need to use the same reference shape and the order of nodes needs to match, e.g.:
```
1---2 "here"
4---3 "there"
InterfaceCell(Line((1,2)), Line((4,3)))
```

# Fields
- `here::AbstractCell`: cell representing the face "here"
- `there::AbstractCell`: cell representing the face "there"
- `nodes`::NTuple: tuple with all node indices in appropriate order: vertex nodes "here", vertex nodes "there", face nodes "here", ...
"""
struct InterfaceCell{shape, Chere, Cthere, N} <: AbstractCell{shape}
    here::Chere
    there::Cthere
    nodes::NTuple{N,Int}

    function InterfaceCell{shape, Chere, Cthere}(here::Chere, there::Cthere) where {shape<:AbstractRefShape, Chere<:AbstractCell, Cthere<:AbstractCell}
        sni = get_sides_and_base_indices(Chere, Cthere)
        nodes = ntuple( i -> sni[i][1] == :here ? here.nodes[sni[i][2]] : there.nodes[sni[i][2]], length(sni))
        new{shape, Chere, Cthere, length(nodes)}(here, there, nodes)
    end
end

function InterfaceCell(here::Chere, there::Cthere) where {Chere<:AbstractCell, Cthere<:AbstractCell}
    @assert getrefshape(here) == getrefshape(there) "For an `InterfaceCell` the underlying cells need to be based on the same shape."
    shape = get_interface_cell_shape(getrefshape(here))
    return InterfaceCell{shape, Chere, Cthere}(here, there)
end

"""
    get_interface_cell_shape(::Type{<:AbstractRefShape})

Return the shape of an interface given a base reference shape.
E.g. given `RefTriangle`, `RefPrism` is returned, meaning two triangles form an interface based on a prism.
"""
get_interface_cell_shape(::Type{RefLine}) = RefQuadrilateral
get_interface_cell_shape(::Type{RefTriangle}) = RefPrism
get_interface_cell_shape(::Type{RefQuadrilateral}) = RefHexahedron


Ferrite.vertices(c::InterfaceCell) = (vertices(c.here)..., vertices(c.there)...)
Ferrite.faces(c::InterfaceCell) = (vertices(c.here), vertices(c.there))
Ferrite.edges(c::InterfaceCell) = (faces(c.here)..., faces(c.there)...)

"""
    get_sides_and_base_indices(c::InterfaceCell)
    get_sides_and_base_indices(::AbstractCell, ::AbstractCell)
    get_sides_and_base_indices(::Type{<:AbstractRefShape}, ::Type{<:AbstractRefShape})

Return a tuple containing tuples of a symbol (:here or :there) and an integer.
The index of the outer tuple represents the node index.
In the inner tuple, the symbol represents the side the node is on 
and the integer represents the nodes index in the base cell.
"""
get_sides_and_base_indices(c::InterfaceCell) = get_sides_and_base_indices(c.here, c.there)
get_sides_and_base_indices(::Chere, ::Cthere) where {Chere <: AbstractCell, Cthere <: AbstractCell} = get_sides_and_base_indices(Chere, Cthere)

get_sides_and_base_indices(::Type{Line}, ::Type{Line}) = ((:here,1), (:here,2), (:there,1), (:there,2))
get_sides_and_base_indices(::Type{QuadraticLine}, ::Type{Line}) = ((:here,1), (:here,2), (:there,1), (:there,2), (:here,3))
get_sides_and_base_indices(::Type{Line}, ::Type{QuadraticLine}) = ((:here,1), (:here,2), (:there,1), (:there,2), (:there,3))
get_sides_and_base_indices(::Type{QuadraticLine}, ::Type{QuadraticLine}) = ((:here,1), (:here,2), (:there,1), (:there,2), (:here,3), (:there,3))

get_sides_and_base_indices(::Type{Triangle}, ::Type{Triangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3))
get_sides_and_base_indices(::Type{QuadraticTriangle}, ::Type{Triangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3), (:here,4), (:here,5), (:here,6))
get_sides_and_base_indices(::Type{Triangle}, ::Type{QuadraticTriangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3), (:there,4), (:there,5), (:there,6))
get_sides_and_base_indices(::Type{QuadraticTriangle}, ::Type{QuadraticTriangle}) = ((:here,1), (:here,2), (:here,3), (:there,1), (:there,2), (:there,3), (:here,4), (:here,5), (:here,6), (:there,4), (:there,5), (:there,6))

get_sides_and_base_indices(::Type{Quadrilateral}, ::Type{Quadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4))
get_sides_and_base_indices(::Type{QuadraticQuadrilateral}, ::Type{Quadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4), (:here,5), (:here,6), (:here,7), (:here,8), (:here,9))
get_sides_and_base_indices(::Type{Quadrilateral}, ::Type{QuadraticQuadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4), (:there,5), (:there,6), (:there,7), (:there,8), (:there,9))
get_sides_and_base_indices(::Type{QuadraticQuadrilateral}, ::Type{QuadraticQuadrilateral}) = ((:here,1), (:here,2), (:here,3), (:here,4), (:there,1), (:there,2), (:there,3), (:there,4), (:here,5), (:here,6), (:here,7), (:here,8), (:there,5), (:there,6), (:there,7), (:there,8), (:here,9), (:there,9))
