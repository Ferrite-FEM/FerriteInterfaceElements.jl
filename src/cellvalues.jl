"""
    InterfaceCellValues([::Type{T},] qr::QuadratureRule, func_ip::InterfaceCellInterpolation, [geom_ip::InterfaceCellInterpolation])

An `InterfaceCellValues` wraps two `CellValues`, one for each face of an `InterfaceCell`.

# Fields
- `ip::InterfaceCellInterpolation`: interpolation on the interface
- `here::CellValues`:  values on face "here"
- `there::CellValues`: values on face "there"
- `basefunctionshere::Vector{Int}`: base function indices on face "here"
- `basefunctionsthere::Vector{Int}`: base function indices on face "there"
"""
struct InterfaceCellValues{IP, CVhere, CVthere} <: AbstractCellValues
    ip::IP
    here::CVhere
    there::CVthere
    basefunctionshere::Vector{Int}
    basefunctionsthere::Vector{Int}

    function InterfaceCellValues(ip::IP, here::CVhere, there::CVthere) where {IP<:Interpolation, CVhere<:CellValues, CVthere<:CellValues}
        @assert here.qr === there.qr "For `InterfaceCellValues` the underlying `CellValues` need to use the same `QuadratureRule`."
        sip = ip isa VectorizedInterpolation ? ip.ip : ip
        basefunctionshere = collect( get_interface_index(sip, :here, i) for i in 1:getnbasefunctions(sip.here))
        basefunctionsthere = collect( get_interface_index(sip, :there, i) for i in 1:getnbasefunctions(sip.there))
        return new{IP, CVhere, CVthere}(ip, here, there, basefunctionshere, basefunctionsthere)
    end
end

function InterfaceCellValues(qr::QuadratureRule,
                             ip::Union{InterfaceCellInterpolation, VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}},
                             gip::Union{InterfaceCellInterpolation, VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}} = default_geometric_interpolation(ip))
    return InterfaceCellValues(Float64, qr, ip, gip)
end

function InterfaceCellValues(::Type{T}, qr::QuadratureRule, ip::InterfaceCellInterpolation) where {T}
    return InterfaceCellValues(T, qr, ip, default_geometric_interpolation(ip))
end
function InterfaceCellValues(::Type{T}, qr::QuadratureRule, ip::VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}) where {T}
    return InterfaceCellValues(T, qr, ip, default_geometric_interpolation(ip.ip))
end

function InterfaceCellValues(::Type{T}, qr::QuadratureRule, 
                             ip::Union{InterfaceCellInterpolation, VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}}, 
                             sgip::InterfaceCellInterpolation{shape}) where {T, sdim, shape <: AbstractRefShape{sdim}}
    return InterfaceCellValues(T, qr, ip, sgip^sdim)
end

function InterfaceCellValues(::Type{T}, qr::QR, ip::IP, gip::VGIP) where {
    T, sdim, rdim, shape <: AbstractRefShape{sdim}, rshape <: AbstractRefShape{rdim},
    QR  <: QuadratureRule{rshape},
    IPhere  <: ScalarInterpolation{rshape},
    IPthere  <: ScalarInterpolation{rshape},
    IP <: InterfaceCellInterpolation{shape, IPhere, IPthere},
    GIPhere <: ScalarInterpolation{rshape},
    GIPthere <: ScalarInterpolation{rshape},
    GIP <: InterfaceCellInterpolation{shape, GIPhere, GIPthere},
    VGIP <: VectorizedInterpolation{sdim, shape, <:Any, GIP}
}
    here = CellValues(T, qr, ip.here, gip.ip.here^sdim)
    there = CellValues(T, qr, ip.there, gip.ip.there^sdim)
    return InterfaceCellValues(ip, here, there)
end

function InterfaceCellValues(::Type{T}, qr::QR, ip::VIP, gip::VGIP) where {
    T, sdim, vdim, rdim, shape <: AbstractRefShape{sdim}, rshape <: AbstractRefShape{rdim},
    QR  <: QuadratureRule{rshape},
    IPhere  <: ScalarInterpolation{rshape},
    IPthere  <: ScalarInterpolation{rshape},
    IP <: InterfaceCellInterpolation{shape, IPhere, IPthere},
    VIP <: VectorizedInterpolation{vdim, shape, <:Any, IP},
    GIPhere <: ScalarInterpolation{rshape},
    GIPthere <: ScalarInterpolation{rshape},
    GIP <: InterfaceCellInterpolation{shape, GIPhere, GIPthere},
    VGIP <: VectorizedInterpolation{sdim, shape, <:Any, GIP}
}
    here = CellValues(T, qr, ip.ip.here^vdim, gip.ip.here^sdim)
    there = CellValues(T, qr, ip.ip.there^vdim, gip.ip.there^sdim)
    return InterfaceCellValues(ip, here, there)
end

Ferrite.reinit!(cv::InterfaceCellValues, cc::CellCache) = reinit!(cv, cc.coords)

function Ferrite.reinit!(cv::InterfaceCellValues, x::AbstractVector{Vec{sdim,T}}) where {sdim, T}
    reinit!(cv.here, @view x[cv.basefunctionshere])
    reinit!(cv.there, @view x[cv.basefunctionsthere])
    return nothing
end

Ferrite.getnbasefunctions(cv::InterfaceCellValues) = getnbasefunctions(cv.here) + getnbasefunctions(cv.there)

Ferrite.getngeobasefunctions(cv::InterfaceCellValues) = getngeobasefunctions(cv.here) + getngeobasefunctions(cv.there)

Ferrite.getnquadpoints(cv::InterfaceCellValues) = getnquadpoints(cv.here)

"""
    get_side_and_baseindex(cv::InterfaceCellValues, i::Integer)

For an `::InterfaceCellValues`: given the base function index `i` return the side (`:here` or `:there`)
the base function belongs to and the corresponding index for the `CellValues` belonging to that side.
"""
function get_side_and_baseindex(cv::InterfaceCellValues, i::Integer)
    ip = cv.ip
    nvhere, nvthere = get_n_dofs_on_side(vertexdof_indices, ip, :here), get_n_dofs_on_side(vertexdof_indices, ip, :there)
    nfhere, nfthere = get_n_dofs_on_side(facedof_interior_indices, ip, :here), get_n_dofs_on_side(facedof_interior_indices, ip, :there)
    nchere, ncthere = get_n_dofs_on_side(celldof_interior_indices, ip, :here), get_n_dofs_on_side(celldof_interior_indices, ip, :there)
    if i ≤ nvhere
        return :here, i
    elseif i ≤ nvhere + nvthere
        return :there, i - nvhere
    elseif i ≤ nvhere + nvthere + nfhere
        return :here, i - nvthere
    elseif i ≤ nvhere + nvthere + nfhere + nfthere
        return :there, i - nvhere - nfhere
    elseif i ≤ nvhere + nvthere + nfhere + nfthere + nchere
        return :here, i - nvthere - nfthere
    elseif i ≤ nvhere + nvthere + nfhere + nfthere + nchere + ncthere
        return :there, i - nvhere - nfhere - nchere
    end
    throw(ArgumentError("No baseindex for interface index $(i) for interpolation $(ip)."))
end

"""
    get_base_value(get_value::Function, cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)

Return a value from an `::InterfaceCellValues` by specifing:
- `get_value`: function specifing which kind of value, e.g. `shape_value`
- `qp`: index of the quadrature point
- `i`: index of the base function
- `here`: side of the interface, where `true` means "here" and `false` means "there".
"""
function get_base_value(get_value::Function, cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)
    side, baseindex = get_side_and_baseindex(cv, i)
    if side == :here && here
        return get_value(cv.here, qp, baseindex)
    elseif side == :there && ! here
        return get_value(cv.there, qp, baseindex)
    else
        return nothing
    end
end

"""
    shape_value(cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)

Return the value of shape function `i` evaluated in quadrature point `qp`
on side `here`, where `true` means "here" and `false` means "there".
"""
function Ferrite.shape_value(cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)
    val = get_base_value(shape_value, cv, qp, i, here)
    if isnothing(val)
        return zero(shape_value_type(cv))
    end
    return val
end

"""
    shape_gradient(cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)

Return the gradient of shape function `i` evaluated in quadrature point `qp`
on side `here`, where `true` means "here" and `false` means "there".
"""
function Ferrite.shape_gradient(cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)
    grad = get_base_value(shape_gradient, cv, qp, i, here)
    if isnothing(grad)
        return zero(shape_gradient_type(cv))
    end
    return grad
end

"""
    shape_value_average(cv::InterfaceCellValues, qp::Int, i::Int)

Return the value of shape function `i` evaluated in quadrature point `qp`
for computing the average value on an interface.
"""
function Ferrite.shape_value_average(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = get_side_and_baseindex(cv, i)
    return side == :here ? shape_value(cv.here, qp, baseindex) / 2 : shape_value(cv.there, qp, baseindex) / 2
end

"""
    shape_gradient_average(cv::InterfaceCellValues, qp::Int, i::Int)

Return the gradient of shape function `i` evaluated in quadrature point `qp`
for computing the average gradient on an interface.
"""
function Ferrite.shape_gradient_average(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = get_side_and_baseindex(cv, i)
    return side == :here ? shape_gradient(cv.here, qp, baseindex) / 2 : shape_gradient(cv.there, qp, baseindex) / 2
end

"""
    shape_value_jump(cv::InterfaceCellValues, qp::Int, i::Int)

Return the value of shape function `i` evaluated in quadrature point `qp`
for computing the value jump on an interface.
"""
function Ferrite.shape_value_jump(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = get_side_and_baseindex(cv, i)
    return side == :here ? -shape_value(cv.here, qp, baseindex) : shape_value(cv.there, qp, baseindex)
end

"""
    shape_gradient_jump(cv::InterfaceCellValues, qp::Int, i::Int)

Return the gradient of shape function `i` evaluated in quadrature point `qp`
for computing the gradient jump on an interface.
"""
function Ferrite.shape_gradient_jump(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = get_side_and_baseindex(cv, i)
    return side == :here ? -shape_gradient(cv.here, qp, baseindex) : shape_gradient(cv.there, qp, baseindex)
end

shape_value_type(cv::InterfaceCellValues) = shape_value_type(cv.here)
shape_gradient_type(cv::InterfaceCellValues) = shape_gradient_type(cv.here)

"""
    function_value(cv::InterfaceCellValues, qp::Int, u::AbstractVector, here::Bool)

Compute the value of the function in a quadrature point on side `here`,
where `true` means "here" and `false` means "there".
`u` is a vector with values for the degrees of freedom.
"""
function Ferrite.function_value(cv::InterfaceCellValues, qp::Int, u::AbstractVector, here::Bool, dof_range = eachindex(u))
    nbf = getnbasefunctions(cv)
    length(dof_range) == nbf || throw_incompatible_dof_length(length(dof_range), nbf)
    @boundscheck checkbounds(u, dof_range)
    @boundscheck checkquadpoint(cv, qp)
    val = function_value_init(cv, u)
    @inbounds for (i, j) in pairs(dof_range)
        val += shape_value(cv, qp, i, here) * u[j]
    end
    return val
end

"""
    function_gradient(cv::InterfaceCellValues, qp::Int, u::AbstractVector, here::Bool)

Compute the gradient of the function in a quadrature point on side `here`,
where `true` means "here" and `false` means "there".
`u` is a vector with values for the degrees of freedom.
"""
function Ferrite.function_gradient(cv::InterfaceCellValues, qp::Int, u::AbstractVector, here::Bool, dof_range = eachindex(u))
    nbf = getnbasefunctions(cv)
    length(dof_range) == nbf || throw_incompatible_dof_length(length(dof_range), nbf)
    @boundscheck checkbounds(u, dof_range)
    @boundscheck checkquadpoint(cv, qp)
    grad = function_gradient_init(cv, u)
    @inbounds for (i, j) in pairs(dof_range)
        grad += shape_gradient(cv, qp, i, here) * u[j]
    end
    return grad
end

"""
    function_value_average(cv::InterfaceCellValues, qp::Int, u::AbstractVector)

Compute the average value of the function in a quadrature point.
"""
function Ferrite.function_value_average(cv::InterfaceCellValues, qp::Int, u::AbstractVector)
    return (function_value(cv, qp, u, true) + function_value(cv, qp, u, false))/2
end

"""
    function_gradient_average(cv::InterfaceCellValues, qp::Int, u::AbstractVector)

Compute the average gradient of the function in a quadrature point.
"""
function Ferrite.function_gradient_average(cv::InterfaceCellValues, qp::Int, u::AbstractVector)
    return (function_gradient(cv, qp, u, true) + function_gradient(cv, qp, u, false)) / 2
end

"""
    getdetJdV_average(cv::InterfaceCellValues, qp::Int)

Return the average of the product between the determinant of the Jacobian on each side of the
interface and the quadrature point weight for the given quadrature point: ``\\det(J(\\mathbf{x})) w_q``.

This value is typically used when integrating a function on the mid-plane of an interface element.
"""
function getdetJdV_average(cv::InterfaceCellValues, qp::Int)
    return (getdetJdV(cv.here, qp) + getdetJdV(cv.there, qp)) / 2
end

"""
    function_value_jump(cv::InterfaceCellValues, qp::Int, u::AbstractVector)

Compute the jump of the function value in a quadrature point.
"""
function Ferrite.function_value_jump(cv::InterfaceCellValues, qp::Int, u::AbstractVector)
    return function_value(cv, qp, u, false) - function_value(cv, qp, u, true)
end

"""
    function_gradient_jump(cv::InterfaceCellValues, qp::Int, u::AbstractVector)

Compute the jump of the function gradient in a quadrature point.
"""
function Ferrite.function_gradient_jump(cv::InterfaceCellValues, qp::Int, u::AbstractVector)
    return function_gradient(cv, qp, u, false) - function_gradient(cv, qp, u, true)
end
