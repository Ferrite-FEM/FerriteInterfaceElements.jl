"""
    InterfaceCell(here::AbstractCell, there::AbstractCell) <: AbstractCell

An `InterfaceCell` is a cell based on two cells of lower dimension representing the two faces.
The two base cells need to be of the same type and the order of nodes needs to match, e.g.:
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
Ferrite.faces(c::InterfaceCell) = (vertices(c.here), vertices(c.there))
Ferrite.edges(c::InterfaceCell) = (faces(c.here)..., faces(c.there)...)

######################################################################
# Inserting cells into a grid
######################################################################

"""
    insert_interfaces(grid, domain_names; topology=ExclusiveTopology(grid))

Return a new grid with `InterfaceCell`s inserted betweenthe domains defined by `domain_names`.
The new grid provides additional cell sets. The set `"interfaces"` contains all new `InterfaceCell`s
and two sets are provided for each combination of domain names: `"domain1-domain2-inertface"` and 
`"domain2-domain1-inertface"` both using the same `Set`.
"""
function insert_interfaces(grid, domain_names; topology=ExclusiveTopology(grid))
    cellsets = Dict(name => getcellset(grid, name) for name in domain_names)
    pairs = Set{Pair{String,Set{Int}}}() # prepare a cellset for each interface between domains
    for i in eachindex(domain_names)
        for j in i+1:length(domain_names)
            set = Set{Int}()
            # add the set with two names to allow for both orders of domain names
            push!(pairs, "$(domain_names[i])-$(domain_names[j])-interface" => set)
            push!(pairs, "$(domain_names[j])-$(domain_names[i])-interface" => set)
        end
    end
    interfacesets = Dict(pairs)
    node_mapping = Dict(name => Dict{Int, Int}() for name in domain_names)

    nodes = copy(grid.nodes)
    cells_generic = Vector{Ferrite.AbstractCell}(grid.cells) # copies

    while !isempty(cellsets)
        name, cellset = pop!(cellsets)
        for cellid in cellset
            cell = getcells(grid, cellid)
            for faceid in 1:nfaces(cell)
                face_neighbors = getneighborhood(topology, grid, FaceIndex(cellid, faceid))
                isempty(face_neighbors) && continue
                (_cellid, _faceid) = only(face_neighbors) # should only ever be one neighboring face
                for (_name, _cellset) in cellsets # only "other" cellsets left
                    if _cellid ∈ _cellset # interface detected
                        facenodeids = Ferrite.faces(cell)[faceid] # original nodeids
                        for nodeid in facenodeids
                            new_nodeid = get(node_mapping[name], nodeid, nothing)
                            _new_nodeid = get(node_mapping[_name], nodeid, nothing)
                            # generate missing duplicate nodes
                            if isnothing(new_nodeid) && isnothing(_new_nodeid)
                                # decide for main side, genererate new node
                                node_mapping[name][nodeid] = nodeid # main node
                                # new node
                                push!(nodes, nodes[nodeid])
                                node_mapping[_name][nodeid] = length(nodes) # main node
                            elseif !isnothing(new_nodeid) && isnothing(_new_nodeid)
                                # node has been duplicated at least once before, so it already has a main
                                # generate a new node for (_name, _cellset)
                                push!(nodes, nodes[nodeid])
                                node_mapping[_name][nodeid] = length(nodes) # main node
                            elseif isnothing(new_nodeid) && !isnothing(_new_nodeid)
                                # node has been duplicated at least once, so it already has a main grain
                                # generate a new node for (name, cellset)
                                push!(nodes, nodes[nodeid])
                                node_mapping[name][nodeid] = length(nodes) # main node
                            end
                        end
                        new_nodeids = Tuple(node_mapping[name][i] for i in facenodeids)
                        _new_nodeids = Tuple(node_mapping[_name][i] for i in facenodeids)
                        # generate new cell
                        @assert typeof(cell) == typeof(getcells(grid, _cellid))
                        interface_cell = create_interface_cell(typeof(cell), new_nodeids, _new_nodeids)
                        push!(cells_generic, interface_cell)
                        push!(interfacesets["$(name)-$(_name)-interface"], length(cells_generic))
                    end
                end
            end
        end
    end

    # better typing of cells vector
    cell_type = Union{(typeof.(unique(typeof,cells_generic))...)}
    cells = convert(Array{cell_type}, cells_generic)

    # adjust original cells to new node numbering
    for name in domain_names
        cellset = getcellset(grid, name)
        for cellid in cellset
            cell = getcells(grid, cellid)
            cells[cellid] = typeof(cell)(map(n->get(node_mapping[name], n, n), cell.nodes))
        end
    end

    new_cellsets = merge(grid.cellsets, Dict("interfaces"=>Set((length(grid.cells)+1):length(cells))), interfacesets)
    # nodesets might no longer be valid
    new_nodesets = Dict{String, Set{Int}}()
    new_grid = Grid(cells, nodes, new_cellsets, new_nodesets, grid.facesets, grid.edgesets, grid.vertexsets, grid.boundary_matrix)
    return new_grid
end

"""
    create_interface_cell(::Type{C}, nodes_here, nodes_there) where {C}

Return a suitable `InterfaceCell` connecting the faces with `nodes_here` and `nodes_there`.
"""
function create_interface_cell(::Type{C}, nodes_here, nodes_there) where {C}
    Cbase = get_interface_base_cell_type(C)
    return InterfaceCell(Cbase(nodes_here), Cbase(nodes_there))
end
"""
    get_interface_base_cell_type(::Type{<:AbstractCell})

Return a suitable base type for connecting two cells of given type with an `InterfaceCell`.
"""
get_interface_base_cell_type(::Type{Triangle}) = Line
get_interface_base_cell_type(::Type{QuadraticTriangle}) = QuadraticLine
get_interface_base_cell_type(::Type{Quadrilateral}) = Line
get_interface_base_cell_type(::Type{QuadraticQuadrilateral}) = QuadraticLine
get_interface_base_cell_type(::Type{Tetrahedron}) = Triangle
get_interface_base_cell_type(::Type{Hexahedron}) = Quadrilateral
