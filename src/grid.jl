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
    relevant_domains, cellsets = _prepare_cellsets(grid, interfaces)
    interfacesets = _prepare_interfacesets(interfaces)
    node_mapping = Dict(name => Dict{Int, Int}() for name in keys(cellsets))

    nodes = copy(grid.nodes)
    cells_generic = Vector{Ferrite.AbstractCell}(grid.cells) # copies

    for (name, domain_h, domain_t) in interfaces
        cellset_h  = cellsets[domain_h]
        cellset_t = cellsets[domain_t]
        for cellid_h in cellset_h
            cell_h = getcells(grid, cellid_h)
            for facetid_h in 1:length(facets(cell_h))
                facet_neighbors = getneighborhood(topology, grid, FacetIndex(cellid_h, facetid_h))
                isempty(facet_neighbors) && continue
                (cellid_t, facetid_t) = only(facet_neighbors) # should only ever be one neighboring face
                if cellid_t in cellset_t # relevant interface detected
                    facetnodeids = Ferrite.facets(cell_h)[facetid_h] # original nodeids
                    for nodeid in facetnodeids
                        new_nodeid_h = get(node_mapping[domain_h], nodeid, nothing)
                        new_nodeid_t = get(node_mapping[domain_t], nodeid, nothing)
                        # generate missing duplicate nodes
                        if isnothing(new_nodeid_h) && isnothing(new_nodeid_t)
                            # decide for main side, genererate new node
                            node_mapping[domain_h][nodeid] = nodeid # main node
                            # new node
                            push!(nodes, nodes[nodeid])
                            node_mapping[domain_t][nodeid] = length(nodes) # main node
                        elseif !isnothing(new_nodeid_h) && isnothing(new_nodeid_t)
                            # node has been duplicated at least once before, so it already has a main
                            # generate a new node for (_name, _cellset)
                            push!(nodes, nodes[nodeid])
                            node_mapping[domain_t][nodeid] = length(nodes) # main node
                        elseif isnothing(new_nodeid_h) && ! isnothing(new_nodeid_t)
                            # node has been duplicated at least once, so it already has a main grain
                            # generate a new node for (name, cellset)
                            push!(nodes, nodes[nodeid])
                            node_mapping[domain_h][nodeid] = length(nodes) # main node
                        end
                    end
                    new_nodeids_h = Tuple(node_mapping[domain_h][i] for i in facetnodeids)
                    new_nodeids_t = Tuple(node_mapping[domain_t][i] for i in facetnodeids)
                    # generate new cell
                    cell_t = getcells(grid, cellid_t)
                    interface_cell = create_interface_cell(typeof(cell_h), typeof(cell_t), new_nodeids_h, new_nodeids_t)
                    push!(cells_generic, interface_cell)
                    _add_interfacecell!(interfacesets, length(cells_generic), name)
                end
            end
        end
    end

    # better typing of cells vector
    cell_type = Union{(typeof.(unique(typeof,cells_generic))...)}
    cells = convert(Array{cell_type}, cells_generic)

    # adjust original cells to new node numbering
    for domain in relevant_domains
        cellset = getcellset(grid, domain)
        for cellid in cellset
            cell = getcells(grid, cellid)
            cells[cellid] = typeof(cell)(map(n -> get(node_mapping[domain], n, n), cell.nodes))
        end
    end

    new_cellsets = merge(grid.cellsets, Dict("interfaces" => OrderedSet((length(grid.cells)+1):length(cells))), interfacesets)
    # nodesets might no longer be valid
    new_nodesets = Dict{String, OrderedSet{Int}}()
    new_grid = Grid(cells, nodes, new_cellsets, new_nodesets, grid.facetsets, grid.vertexsets)
    return new_grid
end

function _prepare_cellsets(grid::Grid, interfaces::Vector{Tuple{T,String,String}}) where {T<:Union{String,Tuple{String,String}}}
    relevant_domains = Set{String}()
    for (_, domain_h, domain_t) in interfaces
        push!(relevant_domains, domain_h)
        push!(relevant_domains, domain_t)
    end
    return relevant_domains, Dict(name => getcellset(grid, name) for name in relevant_domains)
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
