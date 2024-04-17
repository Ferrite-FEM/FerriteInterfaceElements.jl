_ndofs(dofs::Tuple) = sum(idxs -> length(idxs), dofs)
_nvertexdofs(ip) = _ndofs(Ferrite.vertexdof_indices(ip))
_nedgedofs(ip) = _ndofs(Ferrite.edgedof_interior_indices(ip))
_nfacedofs(ip) = _ndofs(Ferrite.facedof_interior_indices(ip))
_ncelldofs(ip) = _ndofs(Ferrite.celldof_interior_indices(ip))

# The constructors of `InterpolationInfo` require a `VectorizedInterpolation{InterfaceCellInterpolation}`
# and not an `InterfaceCellInterpolation{VectorizedInterpolation,VectorizedInterpolation}`.
# To create a `VectorizedInterpolation{InterfaceCellInterpolation}`, `InterfaceCellInterpolation` needs to be a `ScalarInterpolation`.

"""
    InterfaceCellInterpolation(ip::ScalarInterpolation) <: ScalarInterpolation

An `InterfaceCellInterpolation` is based on a regular interpolation which will be applied 
to both faces of an `InterfaceCell`.

# Fields
- `base::ScalarInterpolation`: base interpolation
"""
struct InterfaceCellInterpolation{shape, order, IP} <: ScalarInterpolation{shape, order}
    base::IP

    function InterfaceCellInterpolation(base::IP) where {order, IP<:ScalarInterpolation{<:Any, order}}
        shape = get_interface_cell_shape(getrefshape(base))
        return new{shape, order, IP}(base)
    end
end

Ferrite.getnbasefunctions(ip::InterfaceCellInterpolation) = 2*getnbasefunctions(ip.base)

Ferrite.adjust_dofs_during_distribution(ip::InterfaceCellInterpolation) = Ferrite.adjust_dofs_during_distribution(ip.base)

Ferrite.n_components(ip::InterfaceCellInterpolation) = n_components(ip.base)

Ferrite.getorder(ip::InterfaceCellInterpolation) = getorder(ip.base)
Ferrite.getorder(ip::VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}) = getorder(ip.ip)

function Ferrite.vertexdof_indices(ip::InterfaceCellInterpolation)
    here  = Ferrite.vertexdof_indices(ip.base)
    there = map(v -> map(d -> d + _nvertexdofs(ip.base), v), here)
    return (here..., there...)
end

function Ferrite.edgedof_interior_indices(ip::InterfaceCellInterpolation{<:AbstractRefShape{3}})
    basedofs = Ferrite.facedof_interior_indices(ip.base)
    offset = _nvertexdofs(ip.base)
    here  = map(v -> map(d -> d + offset, v), basedofs)
    offset += _nfacedofs(ip.base)
    there = map(v -> map(d -> d + offset, v), basedofs)
    return (here..., there...)
end

function Ferrite.facedof_indices(ip::InterfaceCellInterpolation{<:AbstractRefShape{dim}}) where {dim}
    basedofs = Ferrite.celldof_indices(ip.base)
    offset = _nvertexdofs(ip.base) + (dim == 3 ? _nedgedofs(ip.base) : 0)
    here  = map(v -> map(d -> d + offset, v), basedofs)
    offset += _ncelldofs(ip.base)
    there = map(v -> map(d -> d + offset, v), basedofs)
    return (here, there)
end

function Ferrite.facedof_interior_indices(ip::InterfaceCellInterpolation{<:AbstractRefShape{dim}}) where {dim}
    basedofs = Ferrite.celldof_interior_indices(ip.base)
    offset = _nvertexdofs(ip.base) + (dim == 3 ? _nedgedofs(ip.base) : 0)
    here  = map(v -> map(d -> d + offset, v), basedofs)
    offset += _ncelldofs(ip.base)
    there = map(v -> map(d -> d + offset, v), basedofs)
    return (here, there)
end

function Ferrite.default_interpolation(::Type{InterfaceCell{shape, C}}) where {shape<:AbstractRefShape, C} 
    return InterfaceCellInterpolation(default_interpolation(C))
end

function Ferrite.default_geometric_interpolation(ip::InterfaceCellInterpolation{<:AbstractRefShape{sdim}, IP}) where {sdim, IP<:ScalarInterpolation}
    return InterfaceCellInterpolation(ip.base)^sdim
end
function Ferrite.default_geometric_interpolation(ip::VectorizedInterpolation{<:Any, <:Any, <:Any, <:InterfaceCellInterpolation})
    return ip
end

"""
    get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)

Return the base function index for an `InterfaceCellInterpolation` given a `side` (`:here` or `:there`)
and the local base function index `i` on that face.
"""
function get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)
    nv = _nvertexdofs(ip.base)
    nf = _nfacedofs(ip.base)
    nc = _ncelldofs(ip.base)
    if side == :here
        if i ≤ nv
            return i
        elseif i ≤ nv + nf
            return i + nv
        elseif i ≤ nv + nf + nc
            return i + nv + nf
        end
        throw(ArgumentError("No interface index for base index $(i) on side $(side) for interpolation $(ip)."))
    elseif side == :there
        if i ≤ nv
            return i + nv
        elseif i ≤ nv + nf
            return i + nv + nf
        elseif i ≤ nv + nf + nc
            return i + nv + nf + nc
        end
        throw(ArgumentError("No interface index for base index $(i) on side $(side) for interpolation $(ip)."))
    end
    throw(ArgumentError("Interface side must be defined by `:here` or `:there`."))
end

"""
    get_interface_dof_indices(get_dofs::Function, ip::InterfaceCellInterpolation) 

Return a tuple of tuples with DOF indices for different entities (vertices, faces, etc.).
The function `get_dofs` specifies which DOFs are considered, e.g. by passing `vertexdof_indices`.
"""
function get_interface_dof_indices(get_dofs::Function, ip::InterfaceCellInterpolation) 
    here  = get_interface_dof_indices(get_dofs, ip, :here)
    there = get_interface_dof_indices(get_dofs, ip, :there)
    return (here..., there...)
end
function get_interface_dof_indices(get_dofs::Function, ip::InterfaceCellInterpolation, side::Symbol)
    basedofs = get_dofs(getproperty(ip, side))
    if isempty(basedofs)
        return (Tuple{}(),)
    else
        return broadcast.(get_interface_index, ((ip,),), ((side,),), basedofs)
    end
end
