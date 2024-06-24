# Cell Values

```@docs
InterfaceCellValues
Ferrite.shape_value(cv::InterfaceCellValues, qp::Int, i::Int; here::Bool)
Ferrite.shape_gradient(cv::InterfaceCellValues, qp::Int, i::Int; here::Bool)
Ferrite.function_value(cv::InterfaceCellValues, qp::Int, u::AbstractVector, dof_range=eachindex(u); here::Bool)
Ferrite.function_gradient(cv::InterfaceCellValues, qp::Int, u::AbstractVector, dof_range=eachindex(u); here::Bool)
shape_value_average
shape_gradient_average
shape_value_jump
shape_gradient_jump
function_value_average
function_gradient_average
function_value_jump
function_gradient_jump
getdetJdV_average
get_side_and_baseindex
```