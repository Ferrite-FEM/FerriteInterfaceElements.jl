# The constructors of `InterpolationInfo` require a `VectorizedInterpolation{InterfaceCellInterpolation}`
# and not an `InterfaceCellInterpolation{VectorizedInterpolation,VectorizedInterpolation}`.
# To create a `VectorizedInterpolation{InterfaceCellInterpolation}`, `InterfaceCellInterpolation` needs to be a `ScalarInterpolation`.

"""
    InterfaceCellInterpolation(here::ScalarInterpolation, there::ScalarInterpolation) <: ScalarInterpolation

An `InterfaceCellInterpolation` is an interpolation based on two interpolations on the faces of an `InterfaceCell`.
If only one interpolation is given, it will be used for both faces.

# Fields
- `here::ScalarInterpolation`: interpolation on the face "here"
- `there::ScalarInterpolation`: interpolation on the face "there"
"""
struct InterfaceCellInterpolation{shape, IPhere, IPthere} <: ScalarInterpolation{shape, Nothing}
    here::IPhere
    there::IPthere

    function InterfaceCellInterpolation(here::IPhere, there::IPthere) where {IPhere<:ScalarInterpolation, IPthere<:ScalarInterpolation}
        @assert getrefshape(here) == getrefshape(there) "For an `InterfaceCellInterpolation` the underlying interpolations need to be based on the same shape."
        shape = get_interface_cell_shape(getrefshape(here))
        return new{shape, IPhere, IPthere}(here, there)
    end
end

function InterfaceCellInterpolation(ip::ScalarInterpolation)
    return InterfaceCellInterpolation(ip, ip)
end

Ferrite.getnbasefunctions(ip::InterfaceCellInterpolation) = getnbasefunctions(ip.here) + getnbasefunctions(ip.there)

"""
    get_n_dofs_on_side(get_dofs::Function, ip::InterfaceCellInterpolation, side::Symbol)

Return the number of DOFs on a `side` (`:here` or `:there`) of an `InterfaceCellInterpolation`.
The function `get_dofs` specifies which DOFs are considered, e.g. by passing `vertexdof_indices`.
"""
function get_n_dofs_on_side(get_dofs::Function, ip::InterfaceCellInterpolation, side::Symbol)
    baseip = getproperty(ip, side)
    return length(get_dofs(baseip)) > 0 ? sum(length.(get_dofs(baseip))) : 0
end
function get_n_dofs_on_side(get_dofs::Function, ip::VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}, side::Symbol)
    baseip = getproperty(ip.ip, side)
    return length(get_dofs(baseip)) > 0 ? sum(length.(get_dofs(baseip)))*n_components(ip) : 0
end

"""
    get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)

Return the base function index for an `InterfaceCellInterpolation` given a `side` (`:here` or `:there`)
and the local base function index `i` on that face.
"""
function get_interface_index(ip::InterfaceCellInterpolation, side::Symbol, i::Integer)
    nvhere, nvthere = get_n_dofs_on_side(vertexdof_indices, ip, :here), get_n_dofs_on_side(vertexdof_indices, ip, :there)
    nfhere, nfthere = get_n_dofs_on_side(facedof_interior_indices, ip, :here), get_n_dofs_on_side(facedof_interior_indices, ip, :there)
    nchere, ncthere = get_n_dofs_on_side(celldof_interior_indices, ip, :here), get_n_dofs_on_side(celldof_interior_indices, ip, :there)
    if side == :here
        if i ≤ nvhere
            return i
        elseif i ≤ nvhere + nfhere
            return i + nvthere
        elseif i ≤ nvhere + nfhere + nchere
            return i + nvthere + nfthere
        end
        throw(ArgumentError("No interface index for base index $(i) on side $(side) for interpolation $(ip)."))
    elseif side == :there
        if i ≤ nvthere
            return i + nvhere
        elseif i ≤ nvthere + nfthere
            return i + nvhere + nfhere
        elseif i ≤ nvthere + nfthere + ncthere
            return i + nvhere + nfhere + nchere
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

Ferrite.vertexdof_indices(ip::InterfaceCellInterpolation) = get_interface_dof_indices(vertexdof_indices, ip)

function Ferrite.facedof_indices(ip::InterfaceCellInterpolation)
    here = (get_interface_dof_indices(vertexdof_indices, ip, :here)...,
            get_interface_dof_indices(facedof_interior_indices, ip, :here)...,
            get_interface_dof_indices(celldof_interior_indices, ip, :here)...)
    there = (get_interface_dof_indices(vertexdof_indices, ip, :there)...,
             get_interface_dof_indices(facedof_interior_indices, ip, :there)...,
             get_interface_dof_indices(celldof_interior_indices, ip, :there)...)
    return Tuple(vcat(collect( [t...] for t in here )...)), Tuple(vcat(collect( [t...] for t in there )...))
end

Ferrite.facedof_interior_indices(ip::InterfaceCellInterpolation) = get_interface_dof_indices(celldof_interior_indices, ip)

Ferrite.edgedof_indices(ip::InterfaceCellInterpolation) = get_interface_dof_indices(facedof_indices, ip)

Ferrite.edgedof_interior_indices(ip::InterfaceCellInterpolation) = get_interface_dof_indices(facedof_interior_indices, ip)

function Ferrite.adjust_dofs_during_distribution(ip::InterfaceCellInterpolation)
    return adjust_dofs_during_distribution(ip.here) || adjust_dofs_during_distribution(ip.there) # TODO: Is this really the way to do it?
end

Ferrite.n_components(ip::InterfaceCellInterpolation) = n_components(ip.here)

function Ferrite.default_interpolation(::Type{InterfaceCell{shape, Chere, Cthere}}) where {shape<:AbstractRefShape, Chere<:AbstractCell, Cthere<:AbstractCell} 
    return InterfaceCellInterpolation(default_interpolation(Chere), default_interpolation(Cthere))
end

function Ferrite.default_geometric_interpolation(ip::InterfaceCellInterpolation{<:AbstractRefShape{sdim}, IPhere, IPthere}
    ) where {sdim, IPhere<:ScalarInterpolation, IPthere<:ScalarInterpolation}
    return InterfaceCellInterpolation(ip.here, ip.there)^sdim
end
function Ferrite.default_geometric_interpolation(ip::InterfaceCellInterpolation{<:AbstractRefShape{sdim}, IPhere, IPthere}
    ) where {sdim, IPhere<:VectorizedInterpolation, IPthere<:VectorizedInterpolation}
    return InterfaceCellInterpolation(ip.here.ip, ip.there.ip)^sdim
end
function Ferrite.default_geometric_interpolation(ip::VectorizedInterpolation{<:Any, <:Any, <:Any, <:InterfaceCellInterpolation})
    return ip
end

Ferrite.getorder(ip::InterfaceCellInterpolation) = getorder(ip.here) == getorder(ip.there) ? getorder(ip.here) : (getorder(ip.here), getorder(ip.there))
Ferrite.getorder(ip::VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}) = getorder(ip.ip)
