function insert_interfaces(grid, domain_names; topology=ExclusiveTopology(grid))
    cellsets = Dict(name=>getcellset(grid, name) for name in domain_names)
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
                        interface_nodes = (new_nodeids..., _new_nodeids...)
                        interface_cell = create_interface_cell(typeof(cell), interface_nodes)
                        push!(cells_generic, interface_cell)
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

    new_cellsets = merge(grid.cellsets, Dict("interfaces"=>Set((length(grid.cells)+1):length(cells))))
    # nodesets might no longer be valid
    new_nodesets = Dict{String, Set{Int}}()
    new_grid = Grid(cells, nodes, new_cellsets, new_nodesets, grid.facesets, grid.edgesets, grid.vertexsets, grid.boundary_matrix)
    return new_grid
end

# TODO: needs to be replaced by FerriteInterfaceElements cells
function create_interface_cell(::Type{C}, nodes) where C
    interface_celltype = get_interface_celltype(C)
    return interface_celltype(nodes)
end
get_interface_celltype(::Type{Triangle}) = CohesiveQuadrilateral
get_interface_celltype(::Type{QuadraticTriangle}) = CohesiveQuadraticQuadrilateral
get_interface_celltype(::Type{Quadrilateral}) = CohesiveQuadrilateral
get_interface_celltype(::Type{QuadraticQuadrilateral}) = CohesiveQuadraticQuadrilateral
get_interface_celltype(::Type{Tetrahedron}) = CohesiveWedge
get_interface_celltype(::Type{Hexahedron}) = CohesiveHexahedron
