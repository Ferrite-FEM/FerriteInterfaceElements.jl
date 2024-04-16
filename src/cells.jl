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

######################################################################
# Inserting cells into a grid
######################################################################

"""
    find_interfaces(grid::Grid, top::ExclusiveTopology, here::Set, there::Set)

Return a `Set{FaceIndex}` with all faces of cells in `here` that are shared with a cell in `there`
and a `Set{Int}` with all node indices on the interface.
"""
function find_interfaces(grid::Grid, top::ExclusiveTopology, here::Set, there::Set)
    interfaces = Set{FaceIndex}()
    nodes = Set{Int}()
    for cellidₕₑᵣₑ in here
        cellₕₑᵣₑ = getcells(grid, cellidₕₑᵣₑ)
        for faceidₕₑᵣₑ in 1:Ferrite.nfaces(cellₕₑᵣₑ)
            faceindexₕₑᵣₑ = FaceIndex(cellidₕₑᵣₑ, faceidₕₑᵣₑ)
            for (cellid, _) in getneighborhood(top, grid, faceindexₕₑᵣₑ)
                if cellid in there
                    push!(interfaces, faceindexₕₑᵣₑ)
                    union!(nodes, Ferrite.faces(cellₕₑᵣₑ)[faceidₕₑᵣₑ])
                end
            end
        end
    end
    return interfaces, nodes
end

"""
    insert_duplicate_nodes!(grid::Grid, targetnodes::Set{Int})

Add duplicates of all nodes specified by `targetnodes` to `grid.nodes`.
Return a `Dict{Int,Int}` as mapping of a given node index to the index of the duplicate.
"""
function insert_duplicate_nodes!(grid::Grid, targetnodes::Set{Int})
    nodes = getnodes(grid)
    noriginals = length(nodes)
    ntargets = length(targetnodes)
    newnodes = typeof(grid.nodes)(undef, ntargets)
    mapping = sizehint!(Dict{Int,Int}(), ntargets)
    for (i, target) in enumerate(targetnodes)
        newid = noriginals + i
        newnodes[i] = deepcopy(nodes[target])
        mapping[target] = newid
    end
    append!(nodes, newnodes)
    return mapping
end

"""
    face_node_indices(::Type{<:AbstractRefShape})

Return a tuple of tuples containing the local node indices of all faces.
"""
face_node_indices(::Type{Triangle}) = Ferrite.facedof_indices(Lagrange{RefTriangle,1}())
face_node_indices(::Type{QuadraticTriangle}) = Ferrite.facedof_indices(Lagrange{RefTriangle,2}())
face_node_indices(::Type{Quadrilateral}) = Ferrite.facedof_indices(Lagrange{RefQuadrilateral,1}())
face_node_indices(::Type{QuadraticQuadrilateral}) = Ferrite.facedof_indices(Lagrange{RefQuadrilateral,2}())
face_node_indices(::Type{Tetrahedron}) = Ferrite.facedof_indices(Lagrange{RefTetrahedron,1}())
face_node_indices(::Type{QuadraticTetrahedron}) = Ferrite.facedof_indices(Lagrange{RefTetrahedron,2}())
face_node_indices(::Type{Hexahedron}) = Ferrite.facedof_indices(Lagrange{RefHexahedron,1}())
face_node_indices(::Type{QuadraticHexahedron}) = Ferrite.facedof_indices(Lagrange{RefHexahedron,2}())

"""
    create_face_cell(::Type{<:Ferrite.AbstractCell}, cell::Ferrite.AbstractCell, face::Int)

Return a new cell for a given face of an existing cell.
"""
function create_face_cell(::Type{CF}, cell::CC, face::Int) where {CF<:Ferrite.AbstractCell, CC<:Ferrite.AbstractCell}
    return CF( getindex.((cell.nodes,), face_node_indices(CC)[face]) )
end
create_face_cell(cell::Triangle, face::Int) = create_face_cell(Line, cell, face)
create_face_cell(cell::QuadraticTriangle, face::Int) = create_face_cell(QuadraticLine, cell, face)
create_face_cell(cell::Quadrilateral, face::Int) = create_face_cell(Line, cell, face)
create_face_cell(cell::QuadraticQuadrilateral, face::Int) = create_face_cell(QuadraticLine, cell, face)
create_face_cell(cell::Tetrahedron, face::Int) = create_face_cell(Triangle, cell, face)
create_face_cell(cell::QuadraticTetrahedron, face::Int) = create_face_cell(QuadraticTetrahedron, cell, face)
create_face_cell(cell::Hexahedron, face::Int) = create_face_cell(Quadrilateral, cell, face)
create_face_cell(cell::QuadraticHexahedron, face::Int) = create_face_cell(QuadraticQuadrilateral, cell, face)

"""
    create_interface_cell(grid::Grid, interface::FaceIndex, nodemapping::Dict{Int,Int})

Creates an `InterfaceCell` at the specified face in the grid.
The `nodemapping` is used to get the node indices of the required duplicate nodes.
"""
function create_interface_cell(grid::Grid, interface::FaceIndex, nodemapping::Dict{Int,Int})
    cellₕₑᵣₑ  = getcells(grid, interface[1])
    facecellₜₕₑᵣₑ  = create_face_cell(cellₕₑᵣₑ,  interface[2])
    facecellₕₑᵣₑ = typeof(facecellₜₕₑᵣₑ)( getindex.((nodemapping,), facecellₜₕₑᵣₑ.nodes) )
    return InterfaceCell(facecellₕₑᵣₑ, facecellₜₕₑᵣₑ)
end

"""
    create_interface_cells(grid::Grid, interfaces::Set{FaceIndex}, nodemapping::Dict{Int,Int})

Creates `InterfaceCell`s at the specified faces in the grid.
The `nodemapping` is used to get the node indices of the required duplicate nodes.
"""
function create_interface_cells(grid::Grid, interfaces::Set{FaceIndex}, nodemapping::Dict{Int,Int})
    return collect( create_interface_cell(grid, faceindex, nodemapping) for faceindex in interfaces)
end

"""
    collect_cells_to_be_adapted(grid::Grid, top::ExclusiveTopology, interfaces::Set{FaceIndex}, here::Set, there::Set)

Return a `Set{Int}` with the indices of all cells for which existing node indices need 
to be replaced by indices of duplicates.
"""
function collect_cells_to_be_adapted(grid::Grid, top::ExclusiveTopology, interfaces::Set{FaceIndex}, here::Set, there::Set)
    targets = Set{Int}(collect(faceindex[1] for faceindex in interfaces))
    for iface in interfaces
        for cellid in getneighborhood(top, grid, CellIndex(iface[1]))
            # Cells in `here` might need to be changed
			if cellid in here
				push!(targets, cellid)
                for faceid in 1:nfaces(getcells(grid, cellid)) # Find possible InterfaceCell between neighbor and other cell
                    for (checkcell, i) in getneighborhood(top, grid, FaceIndex(cellid, faceid))
                        if getcells(grid, checkcell) isa InterfaceCell
                            j = i == 1 ? 2 : 1
                            neighbor = getneighborhood(top, grid, FaceIndex(checkcell, j))
                            if ! (neighbor[1] in there) 
                                push!(targets, checkcell)
                            end
                        end
                    end
                end
			elseif !(cellid in there) || grid.cells[cellid] isa InterfaceCell
                # Only add neighboring cells whith a shared face
                faceneighbor = false
                for faceid in 1:nfaces(getcells(grid, iface[1]))
                    for (checkcell, _) in getneighborhood(top, grid, FaceIndex(iface[1], faceid))
                        faceneighbor = checkcell == cellid
                        if faceneighbor
                            push!(targets, cellid) 
                            break
                        end
                    end
                    faceneighbor ? break : nothing
                end
            end
        end
    end
    return targets
end

"""
    adapt_cell(cell::AbstractCell, nodemapping::Dict{Int,Int})

Return a new cell where existing node indices habe beenreplaced according to `nodemapping`.
"""
function adapt_cell(cell::C, nodemapping::Dict{Int,Int}) where {C<:Ferrite.AbstractCell}
    nodes = collect(n in keys(nodemapping) ? nodemapping[n] : n for n in cell.nodes)
    return C(Tuple(nodes))
end
function adapt_cell(cell::C, nodemapping::Dict{Int,Int}) where {C<:InterfaceCell}
    here = adapt_cell(cell.here, nodemapping)
    there = adapt_cell(cell.there, nodemapping)
    return InterfaceCell(here, there)
end

"""
    adapt_cells!(grid::Grid, top::ExclusiveTopology, interfaces::Set{FaceIndex}, here::Set, there::Set, nodemapping::Dict{Int,Int})

Change existing node indices as needed for cells at the interface.
"""
function adapt_cells!(grid::Grid, top::ExclusiveTopology, interfaces::Set{FaceIndex}, 
        here::Set, there::Set, nodemapping::Dict{Int,Int})
    targets = collect_cells_to_be_adapted(grid, top, interfaces, here, there)
    for i in targets
        grid.cells[i] = adapt_cell(grid.cells[i], nodemapping)
    end
    return grid
end

"""
    create_interface_cells!(grid::Grid, here, there)

Create `InterfaceCell`s between cells in the sets `here` and `there` (as `String`s or `Set{Int}`).
Duplicate nodes are added to the `grid` and existing cells are adapted,
so that cells in `here` and `there` are no longer connected.
Return a `Vector` with `InterfaceCell`s connecting `here` and `there`.

WARNING: Interfaces which are not closed can lead to ambiguities!
"""
function create_interface_cells!(grid::Grid, here::String, there::String)
    here  = getcellset(grid, here)
    there = getcellset(grid, there)
    return create_interface_cells!(grid, here, there)
end

function create_interface_cells!(grid::Grid, here::Set{Int}, there::Set{Int})
    top = ExclusiveTopology(grid)
    interfaces, nodes = find_interfaces(grid, top, here, there)
    mapping = insert_duplicate_nodes!(grid, nodes)
    interfacecells = create_interface_cells(grid, interfaces, mapping)
    adapt_cells!(grid, top, interfaces, here, there, mapping)
    return interfacecells
end
