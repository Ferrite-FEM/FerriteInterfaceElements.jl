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
