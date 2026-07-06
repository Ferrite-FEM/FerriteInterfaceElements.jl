"""
    InterfaceCellValues([::Type{T},] qr::QuadratureRule, func_ip::InterfaceCellInterpolation, [geom_ip::InterfaceCellInterpolation]; use_same_cv=true)

An `InterfaceCellValues` is based on two `CellValues`: one for each facet of an `InterfaceCell`.
Since one can use the same `CellValues` for both sides, be default the same object is used for better performance.
The keyword argument `use_same_cv` can be set to `false` to disable this behavior, if needed.
Furthermore, the rotation matrix of the midplane in the quadrature points can be computed.
To do so, set `include_R=true`.

# Fields
- `ip::InterfaceCellInterpolation`: interpolation on the interface
- `here::CellValues`:  values for facet "here"
- `there::CellValues`:  values for facet "there"
- `base_indices_here::Vector{Int}`: base function indices on facet "here"
- `base_indices_there::Vector{Int}`: base function indices on facet "there"
- `sides_and_baseindices::Tuple`: side and base function for the base `CellValues` for each base function of the `InterfaceCellValues`
- `R::Union{AbstractVector,Nothing}`: rotation matrix of the midplane in the quadrature points
"""
struct InterfaceCellValues{CV,TR,N} <: AbstractCellValues
    here::CV
    there::CV
    base_indices_here::Vector{Int}
    base_indices_there::Vector{Int}
    sides_and_baseindices::NTuple{N,Tuple{Symbol,Int}}
    R::TR # Union{AbstractVector,Nothing} Rotation matrix in quadrature points

    function InterfaceCellValues(ip::IP, here::CV; use_same_cv::Bool, include_R::Val) where {IP<:InterfaceCellInterpolation, CV<:CellValues}
        N = getnbasefunctions(ip)
        sides_and_baseindices = Tuple( get_side_and_baseindex(ip, i) for i in 1:N )
        base_indices_here  = collect( get_interface_index(ip, :here,  i) for i in 1:getnbasefunctions(ip.base) )
        base_indices_there = collect( get_interface_index(ip, :there, i) for i in 1:getnbasefunctions(ip.base) )
        there = use_same_cv ? here : deepcopy(here)
        R = if include_R === Val(false)
            nothing
        else
            T = eltype(here.detJdV)
            Vector{Tensor{2, Ferrite.getrefdim(ip), T}}(undef, getnquadpoints(here))
        end
        return new{CV,typeof(R),N}(here, there, base_indices_here, base_indices_there, sides_and_baseindices, R)
    end
    function InterfaceCellValues(ip::IP, here::CV; use_same_cv, include_R::Val) where {IP<:VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}, CV<:CellValues}
        N = getnbasefunctions(ip)
        sides_and_baseindices = Tuple( get_side_and_baseindex(ip, i) for i in 1:N )
        ip = ip.ip
        base_indices_here  = collect( get_interface_index(ip, :here,  i) for i in 1:getnbasefunctions(ip.base) )
        base_indices_there = collect( get_interface_index(ip, :there, i) for i in 1:getnbasefunctions(ip.base) )
        there = use_same_cv ? here : deepcopy(here)
        R = if include_R === Val(false)
            nothing
        else
            T = eltype(here.detJdV)
            Vector{Tensor{2, Ferrite.getrefdim(ip), T}}(undef, getnquadpoints(here))
        end
        return new{CV,typeof(R),N}(here, there, base_indices_here, base_indices_there, sides_and_baseindices, R)
    end
end

InterfaceCellValues(qr::QuadratureRule, args...; kwargs...) = InterfaceCellValues(Float64, qr, args...; kwargs...)

function InterfaceCellValues(::Type{T}, qr::QuadratureRule, 
            ip::InterfaceCellInterpolation, 
            ip_geo::VectorizedInterpolation{sdim,<:Any,<:Any,<:InterfaceCellInterpolation} = default_geometric_interpolation(ip); 
            use_same_cv=true, include_R=Val(false), kwargs...) where {T, sdim}
    cv = CellValues(T, qr, ip.base, VectorizedInterpolation{sdim}(ip_geo.ip.base); kwargs...)
    return InterfaceCellValues(ip, cv; use_same_cv=use_same_cv, include_R = isa(include_R, Bool) ? Val(include_R) : include_R)
end
function InterfaceCellValues(::Type{T}, qr::QuadratureRule, 
            ip::VectorizedInterpolation{vdim,<:Any,<:Any,<:InterfaceCellInterpolation}, 
            ip_geo::VectorizedInterpolation{sdim,<:Any,<:Any,<:InterfaceCellInterpolation} = default_geometric_interpolation(ip); 
            use_same_cv=true, include_R=Val(false), kwargs...) where {T, vdim, sdim}
    cv = CellValues(T, qr, VectorizedInterpolation{vdim}(ip.ip.base), 
                           VectorizedInterpolation{sdim}(ip_geo.ip.base); kwargs...)
    return InterfaceCellValues(ip, cv; use_same_cv=use_same_cv, include_R = isa(include_R, Bool) ? Val(include_R) : include_R)
end

function InterfaceCellValues(::Type{T}, qr::QuadratureRule, 
            ip::InterfaceCellInterpolation, 
            ip_geo::InterfaceCellInterpolation{shape}; kwargs...) where {T, dim, shape <: AbstractRefShape{dim}}
    return InterfaceCellValues(T, qr, ip, VectorizedInterpolation{dim}(ip); kwargs...)
end
function InterfaceCellValues(::Type{T}, qr::QuadratureRule, 
            ip::VectorizedInterpolation{<:Any,<:Any,<:Any,<:InterfaceCellInterpolation}, 
            ip_geo::InterfaceCellInterpolation{shape}; kwargs...) where {T, dim, shape <: AbstractRefShape{dim}}
    return InterfaceCellValues(T, qr, ip, VectorizedInterpolation{dim}(ip.ip); kwargs...)
end

function Ferrite.default_geometric_interpolation(::Type{InterfaceCell{shape, C, N}}) where {shape<:AbstractRefShape, C, N} 
    return InterfaceCellInterpolation(geometric_interpolation(C))
end

function Ferrite.default_geometric_interpolation(ip::InterfaceCellInterpolation{<:AbstractRefShape{sdim}, <:Any, IP}) where {sdim, IP<:ScalarInterpolation}
    return VectorizedInterpolation{sdim}(InterfaceCellInterpolation(ip.base))
end
function Ferrite.default_geometric_interpolation(ip::VectorizedInterpolation{<:Any, <:Any, <:Any, <:InterfaceCellInterpolation})
    return ip
end

Ferrite.getnbasefunctions(cv::InterfaceCellValues) = 2*getnbasefunctions(cv.here)

Ferrite.getngeobasefunctions(cv::InterfaceCellValues) = 2*getngeobasefunctions(cv.here)

Ferrite.getnquadpoints(cv::InterfaceCellValues) = getnquadpoints(cv.here)

Ferrite.shape_value_type(cv::InterfaceCellValues) = shape_value_type(cv.here)

Ferrite.shape_gradient_type(cv::InterfaceCellValues) = shape_gradient_type(cv.here)

Ferrite.reinit!(cv::InterfaceCellValues, cc::CellCache) = reinit!(cv, cc.coords)

function Ferrite.reinit!(cv::InterfaceCellValues{CV}, x::AbstractVector{Vec{sdim,T}}) where {sdim, T, CV}
    n_coords_per_side = length(x) ÷ 2
    x_here  = @view x[cv.base_indices_here[1:n_coords_per_side]]
    reinit!(cv.here, x_here)

    if ! (cv.here === cv.there)
        x_there = @view x[cv.base_indices_there[1:n_coords_per_side]]
        reinit!(cv.there, x_there)
    end

    if cv.R !== nothing
        x_there = @view x[cv.base_indices_there]
        for qp in 1:getnquadpoints(cv.here)
            mapping_here  = Ferrite.calculate_mapping(cv.here.geo_mapping,  qp, x_here)
            mapping_there = Ferrite.calculate_mapping(cv.there.geo_mapping, qp, x_there)
            J_here  = Ferrite.getjacobian(mapping_here)
            J_there = Ferrite.getjacobian(mapping_there)
            J = (J_here + J_there) / 2
            cv.R[qp] = _get_R_from_J(J)
        end
    end
    return nothing
end
function _get_R_from_J(J::MixedTensor2{2,1,T}) where T 
    v1 = J[:, 1]
    v2 = Vec{2,T}((-v1[2], v1[1]))
    return Tensor{2,2,T}(((v1/norm(v1))..., (v2/norm(v2))...))
end
function _get_R_from_J(J::MixedTensor2{3,2,T}) where T 
    v1 = J[:, 1]
    v3 = v1 × J[:, 2]
    v2 = v3 × v1
    return Tensor{2,3,T}(((v1/norm(v1))..., (v2/norm(v2))..., (v3/norm(v3))...))
end

get_side_and_baseindex(cv::InterfaceCellValues, i::Integer) = cv.sides_and_baseindices[i]

"""
    get_base_value(get_value::Function, cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)

Return a value from an `::InterfaceCellValues` by specifing:
- `get_value`: function specifing which kind of value, e.g. `shape_value`
- `qp`: index of the quadrature point
- `i`: index of the base function
- `here`: side of the interface, where `true` means "here" and `false` means "there".
"""
function get_base_value(get_value::Function, cv::InterfaceCellValues, qp::Int, i::Int, here::Bool)
    side, baseindex = cv.sides_and_baseindices[i]
    if side == :here && here
        return get_value(cv.here, qp, baseindex)
    elseif side == :there && ! here
        return get_value(cv.there, qp, baseindex)
    end
    return nothing
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
    side, baseindex = cv.sides_and_baseindices[i]
    return side==:here ? shape_value(cv.here, qp, baseindex) / 2 : shape_value(cv.there, qp, baseindex) / 2
end

"""
    shape_gradient_average(cv::InterfaceCellValues, qp::Int, i::Int)

Return the gradient of shape function `i` evaluated in quadrature point `qp`
for computing the average gradient on an interface.
"""
function Ferrite.shape_gradient_average(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = cv.sides_and_baseindices[i]
    return side == :here ? shape_gradient(cv.here, qp, baseindex) / 2 : shape_gradient(cv.there, qp, baseindex) / 2
end

"""
    shape_value_jump(cv::InterfaceCellValues, qp::Int, i::Int)

Return the value of shape function `i` evaluated in quadrature point `qp`
for computing the value jump on an interface.
"""
function Ferrite.shape_value_jump(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = cv.sides_and_baseindices[i]
    return side == :here ? -shape_value(cv.here, qp, baseindex) : shape_value(cv.there, qp, baseindex)
end

"""
    shape_gradient_jump(cv::InterfaceCellValues, qp::Int, i::Int)

Return the gradient of shape function `i` evaluated in quadrature point `qp`
for computing the gradient jump on an interface.
"""
function Ferrite.shape_gradient_jump(cv::InterfaceCellValues, qp::Int, i::Int)
    side, baseindex = cv.sides_and_baseindices[i]
    return side == :here ? -shape_gradient(cv.here, qp, baseindex) : shape_gradient(cv.there, qp, baseindex)
end

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

"""
    midplane_rotation(cv::InterfaceCellValues, qp::Int)

Compute the rotation matrix of the midplane in quadrature point.
"""
function midplane_rotation(cv::InterfaceCellValues, qp::Int)
    return cv.R[qp]
end

function Base.show(io::IO, d::MIME"text/plain", cv::InterfaceCellValues)
    if cv.here === cv.there
        print(io, "InterfaceCellValues based on a single CellValues:\n"); show(io, d, cv.here)
    else
        print(io, "InterfaceCellValues based on two different CellValues: ")
        print(io, "\nHere:\n"); show(io, d, cv.here)
        print(io, "\nThere:\n"); show(io, d, cv.there)
    end
    print(io, "\n+ Additional info:")
    print(io, "\nBase functions here:  ", cv.base_indices_here)
    print(io, "\nBase functions there: ", cv.base_indices_there)
    print(io, "\nSides and indices for base cell values: ", cv.sides_and_baseindices)
end
