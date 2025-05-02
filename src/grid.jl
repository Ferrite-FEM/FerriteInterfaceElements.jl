######################################################################
# Inserting cells into a grid
######################################################################

"""
    insert_interfaces(grid, domain_names; topology=ExclusiveTopology(grid))

Return a new grid with `InterfaceCell`s inserted betweenthe domains defined by `domain_names`.
The new grid provides additional cell sets. The set `"interfaces"` contains all new `InterfaceCell`s
and two sets are provided for each combination of domain names: `"domain1-domain2-interface"` and 
`"domain2-domain1-interface"` both using the same `Set`.
"""
function insert_interfaces(grid, domain_names; topology=ExclusiveTopology(grid))
    cellsets = Dict(name => getcellset(grid, name) for name in domain_names)
    pairs = OrderedSet{Pair{String,OrderedSet{Int}}}() # prepare a cellset for each interface between domains
    for i in eachindex(domain_names)
        for j in i+1:length(domain_names)
            set = OrderedSet{Int}()
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
            for facetid in 1:length(facets(cell))
                facet_neighbors = getneighborhood(topology, grid, FacetIndex(cellid, facetid))
                isempty(facet_neighbors) && continue
                (_cellid, _facetid) = only(facet_neighbors) # should only ever be one neighboring face
                for (_name, _cellset) in cellsets # only "other" cellsets left
                    if _cellid ∈ _cellset # interface detected
                        facetnodeids = Ferrite.facets(cell)[facetid] # original nodeids
                        for nodeid in facetnodeids
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
                        new_nodeids = Tuple(node_mapping[name][i] for i in facetnodeids)
                        _new_nodeids = Tuple(node_mapping[_name][i] for i in facetnodeids)
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

    new_cellsets = merge(grid.cellsets, Dict("interfaces" => OrderedSet((length(grid.cells)+1):length(cells))), interfacesets)
    # nodesets might no longer be valid
    new_nodesets = Dict{String, OrderedSet{Int}}()
    new_grid = Grid(cells, nodes, new_cellsets, new_nodesets, grid.facetsets, grid.vertexsets)
    return new_grid
end

"""
    create_interface_cell(::Type{C}, nodes_here, nodes_there) where {C}

Return a suitable `InterfaceCell` connecting the facets with `nodes_here` and `nodes_there`.
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
