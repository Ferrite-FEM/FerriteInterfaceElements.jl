#Ferrite.facedof_interior_indices(::Lagrange{RefLine}) = Tuple{}() # Hot fix -> TODO: remove as soon as its in Ferrite

# The constructors of `InterpolationInfo` require a `VectorizedInterpolation{InterfaceCellInterpolation}`
# and not an `InterfaceCellInterpolation{VectorizedInterpolation,VectorizedInterpolation}`.
# To create a `VectorizedInterpolation{InterfaceCellInterpolation}`, `InterfaceCellInterpolation` needs to be a `ScalarInterpolation`.

"""
    InterfaceCellInterpolation(ip::ScalarInterpolation) <: ScalarInterpolation

An `InterfaceCellInterpolation` is based on a regular interpolation which will be applied 
to both facets of an `InterfaceCell`.

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

function Ferrite.vertexdof_indices(ip::InterfaceCellInterpolation)
    here  = Ferrite.vertexdof_indices(ip.base)
    offset = _nvertexdofs(ip.base)
    there = map(v -> map(d -> d + offset, v), here)
    return (here..., there...)
end

function Ferrite.edgedof_interior_indices(ip::InterfaceCellInterpolation{<:AbstractRefShape{dim}}) where {dim}
    basedofs = Ferrite.edgedof_interior_indices(ip.base)
    offset = _nvertexdofs(ip.base)
    here  = map(v -> map(d -> d + offset, v), basedofs)
    offset += length(basedofs)
    there = map(v -> map(d -> d + offset, v), basedofs)
    return (here..., there...)
end

function Ferrite.facedof_interior_indices(ip::InterfaceCellInterpolation{<:AbstractRefShape{3}})
    basedofs = Ferrite.facedof_interior_indices(ip.base)
    @assert length(basedofs) == 1 "In 3D an `InterfaceCellInterpolation` is based on a 2D interpolation (with only a single face)."
    if isempty(basedofs[1])
        return (basedofs[1], basedofs[1])
    end
    offset = _nvertexdofs(ip.base) + _nedgedofs(ip.base)
    here  = map(dofs -> map(d -> d + offset, dofs), basedofs)
    offset += length(basedofs)
    there = map(dofs -> map(d -> d + offset, dofs), basedofs)
    return (here, there)
end

#########################################################
# Functions for setting up the cell values
#########################################################

_ndofs(dofs::Tuple) = isempty(dofs) ? 0 : sum(idxs -> length(idxs), dofs)
_nvertexdofs(ip) = _ndofs(Ferrite.vertexdof_indices(ip))
_nedgedofs(ip) = _ndofs(Ferrite.edgedof_interior_indices(ip))
_nfacedofs(ip) = _ndofs(Ferrite.facedof_interior_indices(ip))

_nvertexdofs_perside(ip::InterfaceCellInterpolation) = _nvertexdofs(ip.base)
_nedgedofs_perside(ip::InterfaceCellInterpolation) = _nedgedofs(ip.base)
_nfacedofs_perside(ip::InterfaceCellInterpolation) = _nfacedofs(ip.base)

_nvertexdofs_perside(ip::VectorizedInterpolation{dim,<:Any,<:Any,<:InterfaceCellInterpolation}) where {dim} = dim*_nvertexdofs(ip.ip.base)
_nedgedofs_perside(ip::VectorizedInterpolation{dim,<:Any,<:Any,<:InterfaceCellInterpolation}) where {dim} = dim*_nedgedofs(ip.ip.base)
_nfacedofs_perside(ip::VectorizedInterpolation{dim,<:Any,<:Any,<:InterfaceCellInterpolation}) where {dim}  = dim*_nfacedofs(ip.ip.base)

"""
    get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)

Return the base function index for an `InterfaceCellInterpolation` given a `side` (`:here` or `:there`)
and the local base function index `i` on that face.
"""
function get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)
    nv = _nvertexdofs(ip.base)
    ne = _nedgedofs(ip.base)
    nf = _nfacedofs(ip.base)
    if side == :here
        if i ≤ nv
            return i
        elseif i ≤ nv + ne
            return i + nv
        elseif i ≤ nv + ne + nf
            return i + nv + ne
        end
        throw(ArgumentError("No interface index for base index $(i) on side $(side) for interpolation $(ip)."))
    elseif side == :there
        if i ≤ nv
            return i + nv
        elseif i ≤ nv + ne
            return i + nv + ne
        elseif i ≤ nv + ne + nf
            return i + nv + ne + nf
        end
        throw(ArgumentError("No interface index for base index $(i) on side $(side) for interpolation $(ip)."))
    end
    throw(ArgumentError("Interface side must be defined by `:here` or `:there`."))
end

"""
    get_side_and_baseindex(cv::InterfaceCellInterpolation, i::Integer)
    get_side_and_baseindex(cv::InterfaceCellValues, i::Integer)

For an `InterfaceCellInterpolation`: given the base function index `i` return the side (`:here` or `:there`)
 and the base function index locally on that side.
"""
function get_side_and_baseindex(ip::Union{InterfaceCellInterpolation,
                                          VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}},
                                i::Integer)
    nv = _nvertexdofs_perside(ip)
    ne = _nedgedofs_perside(ip)
    nf = _nfacedofs_perside(ip)
    if i ≤ nv
        return :here, i
    elseif i ≤ 2*nv
        return :there, i - nv
    elseif i ≤ 2*nv + ne
        return :here, i - nv
    elseif i ≤ 2*nv + 2*ne
        return :there, i - nv - ne
    elseif i ≤ 2*nv + 2*ne + nf
        return :here, i - nv - ne
    elseif i ≤ 2*nv + 2*ne + 2*nf
        return :there, i - nv - ne - nf
    end
    throw(ArgumentError("Index $(i) exeeds number of basefunctions."))
end
