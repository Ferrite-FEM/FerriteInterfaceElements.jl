Ferrite.geometric_interpolation(T::Type{InterfaceCell{shape, C, N}}) where {shape<:AbstractRefShape, C, N}  = Ferrite.default_geometric_interpolation(T) # Needed for Ferrite.apply_analytical!

"""
    apply_analytical_to_bulk!(a::AbstractVector, dh::DofHandler{<:Grid{C}},
                      fieldname::Symbol, f::Function, cellset=1:getncells(get_grid(dh))) where {C}

This method catches whether the underlying grid contains `InterfaceCell`s or not.
Then, it calls `Ferrite.apply_analytical!` excluding `InterfaceCell`s.
"""
function apply_analytical_to_bulk!(a::AbstractVector, dh::DofHandler{dim,<:Grid{dim,C}}, # Catch cases possibly involving interface cells
                                   fieldname::Symbol, f::Function, cellset=1:getncells(get_grid(dh))) where {dim,C}
    if isconcretetype(C)
        return Ferrite.apply_analytical!(a, dh, fieldname, f, cellset)
    end
    cellset = OrderedSet{Int}(cellset)
    filter!(i -> ! (dh.grid.cells[i] isa InterfaceCell), cellset)
    return Ferrite.apply_analytical!(a, dh, fieldname, f, cellset)
end
