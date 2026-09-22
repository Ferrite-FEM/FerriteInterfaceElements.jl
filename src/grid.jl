######################################################################
# Inserting cells into a grid
######################################################################
"""
    insert_interfaces(grid::Grid, domain_names::Vector{String}; kwargs...)    
    insert_interfaces(grid::Grid, interfaces::Dict{String,Tuple{String,String}}; kwargs...)
    
Return a new grid with `InterfaceCell`s inserted according to `domain_names` or `interfaces`.
The new grid provides additional cell sets. The set `"interfaces"` contains all new `InterfaceCell`s.

When using `domain_names`, for each combination of the corresponding cellsests interfaces will be inserted
and a new cellset will be provided, which can be accessed using both of the following names:
`"domain1-domain2-interface"` and `"domain2-domain1-interface"`.

When using `interfaces`, each value of the `Dict` defines a pair of cellsets between which interfaces
will be inserted and collected in a new cellset which can be accessed by the corresponding key.

Nodes are only duplicated between cells which are separated by an inserted interface. Where three or
more domains meet in a single node, domains which are not separated by an interface keep sharing the node.
"""
function insert_interfaces(grid::Grid, interfaces::Dict{String,Tuple{String,String}}; kwargs...)
    interfaces = [(name, domains[1], domains[2]) for (name, domains) in pairs(interfaces)]
    return _insert_interfaces(grid, interfaces; kwargs...)
end
function insert_interfaces(grid::Grid, domain_names::Vector{String}; kwargs...)
    interfaces = [begin
        names = ("$(domain_names[i])-$(domain_names[j])-interface", "$(domain_names[j])-$(domain_names[i])-interface")
        (names, domain_names[i], domain_names[j])
        end for i in 1:length(domain_names) for j in i+1:length(domain_names)]
    return _insert_interfaces(grid, interfaces; kwargs...)
end



function _insert_interfaces(grid::Grid, interfaces::Vector{Tuple{T,String,String}}; topology=ExclusiveTopology(grid)) where {T<:Union{String,Tuple{String,String}}}
    cellsets = _prepare_cellsets(grid, interfaces)
    interfacesets = _prepare_interfacesets(interfaces)

    # Collect the facets to be split, in the order in which the interface cells are created
    interfacefacets = Tuple{T, Tuple{Int, Int}, Tuple{Int, Int}}[] # (name, (cellid, facetid) here, (cellid, facetid) there)
    splitfacets = Set{Tuple{Int, Int}}()
    for (name, domain_h, domain_t) in interfaces
        cellset_t = cellsets[domain_t]
        for cellid_h in cellsets[domain_h]
            cell_h = getcells(grid, cellid_h)
            for facetid_h in 1:nfacets(cell_h)
                facet_neighbors = getneighborhood(topology, grid, FacetIndex(cellid_h, facetid_h))
                isempty(facet_neighbors) && continue
                (cellid_t, facetid_t) = only(facet_neighbors) # should only ever be one neighboring facet
                cellid_t in cellset_t || continue
                facet_h, facet_t = (cellid_h, facetid_h), (cellid_t, facetid_t)
                push!(interfacefacets, (name, facet_h, facet_t))
                push!(splitfacets, facet_h, facet_t)
            end
        end
    end

    # Collect the cells around every node on a split facet
    node_cells = Dict{Int, Vector{Int}}()
    for (_, facet_h, _) in interfacefacets
        for nodeid in _facet_node_ids(grid, facet_h)
            get!(node_cells, nodeid, Int[])
        end
    end
    for (cellid, cell) in enumerate(getcells(grid))
        for nodeid in Ferrite.get_node_ids(cell)
            cells = get(node_cells, nodeid, nothing)
            isnothing(cells) || push!(cells, cellid)
        end
    end

    # Group the cells around each split node into components: two cells belong to the same component
    # if they are connected through facets which are not split. Each component gets its own copy of
    # the node. This matters at junctions where three or more domains meet in a single node but
    # interfaces are only requested between some of them: the domains which are not separated by an
    # interface must keep sharing the node.
    node_components = Dict{Int, Dict{Int, Int}}() # nodeid => (cellid => component id)
    for (nodeid, cellids) in node_cells
        parent = Dict(cellid => cellid for cellid in cellids)
        for cellid in cellids
            cell = getcells(grid, cellid)
            for facetid in 1:nfacets(cell)
                (cellid, facetid) in splitfacets && continue
                nodeid in _facet_node_ids(grid, (cellid, facetid)) || continue
                for neighbor in getneighborhood(topology, grid, FacetIndex(cellid, facetid))
                    _union!(parent, cellid, neighbor[1])
                end
            end
        end
        node_components[nodeid] = Dict(cellid => _find!(parent, cellid) for cellid in cellids)
    end

    # Assign node ids to the components: the first component keeps the original node,
    # every other component gets a duplicate.
    nodes = copy(getnodes(grid))
    new_nodeids = Dict{Tuple{Int, Int}, Int}() # (nodeid, component id) => new nodeid
    kept = Set{Int}() # nodes for which a component keeps the original node id
    cells_generic = Vector{Ferrite.AbstractCell}(getcells(grid)) # copies
    for (name, facet_h, facet_t) in interfacefacets
        cellid_h, cellid_t = facet_h[1], facet_t[1]
        facetnodeids = _facet_node_ids(grid, facet_h)
        for nodeid in facetnodeids
            _assign_nodeid!(new_nodeids, kept, nodes, nodeid, node_components[nodeid][cellid_h])
            _assign_nodeid!(new_nodeids, kept, nodes, nodeid, node_components[nodeid][cellid_t])
        end
        nodeids_h = map(n -> new_nodeids[(n, node_components[n][cellid_h])], facetnodeids)
        nodeids_t = map(n -> new_nodeids[(n, node_components[n][cellid_t])], facetnodeids)
        # generate new cell
        interface_cell = create_interface_cell(typeof(getcells(grid, cellid_h)), typeof(getcells(grid, cellid_t)), nodeids_h, nodeids_t)
        push!(cells_generic, interface_cell)
        _add_interfacecell!(interfacesets, length(cells_generic), name)
    end

    # better typing of cells vector
    cell_type = Union{(typeof.(unique(typeof,cells_generic))...)}
    cells = convert(Array{cell_type}, cells_generic)

    # adjust original cells to new node numbering
    for cellid in Set{Int}(Iterators.flatten(values(node_cells)))
        cell = getcells(grid, cellid)
        cells[cellid] = typeof(cell)(map(Ferrite.get_node_ids(cell)) do n
            components = get(node_components, n, nothing)
            return isnothing(components) ? n : new_nodeids[(n, components[cellid])]
        end)
    end

    new_cellsets = merge(Ferrite.getcellsets(grid), Dict("interfaces" => OrderedSet((getncells(grid)+1):length(cells))), interfacesets)
    # nodesets might no longer be valid
    new_nodesets = Dict{String, OrderedSet{Int}}()
    new_grid = Grid(cells, nodes; cellsets=new_cellsets, nodesets=new_nodesets,
                    facetsets=Ferrite.getfacetsets(grid), vertexsets=Ferrite.getvertexsets(grid))
    return new_grid
end

# Node ids of a facet, including any non-vertex nodes (e.g. mid nodes of quadratic cells)
function _facet_node_ids(grid::Grid, (cellid, facetid)::Tuple{Int, Int})
    cell = getcells(grid, cellid)
    facetdofs = Ferrite.facetdof_indices(geometric_interpolation(cell))[facetid]
    return map(i -> Ferrite.get_node_ids(cell)[i], facetdofs)
end

# Union-find over cell ids
function _find!(parent::Dict{Int, Int}, i::Int)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end
function _union!(parent::Dict{Int, Int}, i::Int, j::Int)
    ri, rj = _find!(parent, i), _find!(parent, j)
    ri == rj || (parent[max(ri, rj)] = min(ri, rj))
    return nothing
end

# Return the node id for `nodeid` in the given component, creating a duplicate node if needed.
# The first component that asks for a node keeps the original one.
function _assign_nodeid!(new_nodeids::Dict{Tuple{Int, Int}, Int}, kept::Set{Int}, nodes::Vector, nodeid::Int, component::Int)
    return get!(new_nodeids, (nodeid, component)) do
        if nodeid in kept
            push!(nodes, nodes[nodeid])
            return length(nodes)
        else
            push!(kept, nodeid)
            return nodeid
        end
    end
end

function _prepare_cellsets(grid::Grid, interfaces::Vector{Tuple{T,String,String}}) where {T<:Union{String,Tuple{String,String}}}
    relevant_domains = Set{String}()
    for (_, domain_h, domain_t) in interfaces
        push!(relevant_domains, domain_h)
        push!(relevant_domains, domain_t)
    end
    return Dict(name => getcellset(grid, name) for name in relevant_domains)
end

function _prepare_interfacesets(interfaces::Vector{Tuple{String,String,String}})
    return Dict([ name => OrderedSet{Int}() for (name,_,_) in interfaces ])
end
function _prepare_interfacesets(interfaces::Vector{Tuple{Tuple{String,String},String,String}})
    sets = [OrderedSet{Int}() for _ in interfaces]
    return Dict([ name => set for (set, (names,_,_)) in zip(sets, interfaces) for name in names ])
end

_add_interfacecell!(interfacesets::Dict{String,OrderedSet{Int}}, cellid::Int, name::String) = push!(interfacesets[name], cellid)
_add_interfacecell!(interfacesets::Dict{String,OrderedSet{Int}}, cellid::Int, name::Tuple{String,String}) = push!(interfacesets[name[1]], cellid)


"""
    create_interface_cell(::Type{C₁}, ::Type{C₂}, nodes_h, nodes_t) where {C₁,C₂}

Return a suitable `InterfaceCell` connecting the facets with `nodes_h` and `nodes_t`.
"""
function create_interface_cell(::Type{C₁}, ::Type{C₂}, nodes_h, nodes_t) where {C₁,C₂}
    Cbase = get_interface_base_cell_type(C₁, C₂, nodes_h)
    return InterfaceCell(Cbase(nodes_h), Cbase(nodes_t))
end


"""
    get_interface_base_cell_type(::Type{<:AbstractCell}, ::Type{<:AbstractCell})

Return a suitable base type for connecting two cells of given type with an `InterfaceCell`.
"""
get_interface_base_cell_type(::Type{Triangle}, ::Type{Triangle}, ::Any) = Line
get_interface_base_cell_type(::Type{QuadraticTriangle}, ::Type{QuadraticTriangle}, ::Any) = QuadraticLine
get_interface_base_cell_type(::Type{Quadrilateral}, ::Type{Quadrilateral}, ::Any) = Line
get_interface_base_cell_type(::Type{QuadraticQuadrilateral}, ::Type{QuadraticQuadrilateral}, ::Any) = QuadraticLine
get_interface_base_cell_type(::Type{Tetrahedron}, ::Type{Tetrahedron}, ::Any) = Triangle
get_interface_base_cell_type(::Type{Hexahedron}, ::Type{Hexahedron}, ::Any) = Quadrilateral

get_interface_base_cell_type(::Type{Triangle}, ::Type{Quadrilateral}, ::Any) = Line
get_interface_base_cell_type(::Type{Quadrilateral}, ::Type{Triangle}, ::Any) = Line
get_interface_base_cell_type(::Type{QuadraticTriangle}, ::Type{QuadraticQuadrilateral}, ::Any) = QuadraticLine
get_interface_base_cell_type(::Type{QuadraticQuadrilateral}, ::Type{QuadraticTriangle}, ::Any) = QuadraticLine

function get_interface_base_cell_type(::Type{C₁}, ::Type{C₂}, nodes) where {C₁<:Union{Pyramid,Wedge}, C₂<:Union{Pyramid,Wedge}}
    if length(nodes) == 4
        return Quadrilateral
    elseif length(nodes) == 3
        return Triangle
    end
    throw(ErrorException("No feasible base type for InterfaceCell!"))
    return nothing    
end
