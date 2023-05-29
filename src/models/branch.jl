function branch!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r_from::T,
    voltage_i_from::T,
    voltage_r_to::T,
    voltage_i_to::T,
    current_r_from::AbstractArray{T},
    current_i_from::AbstractArray{T},
    current_r_to::AbstractArray{T},
    current_i_to::AbstractArray{T},
    branch::BranchWrapper,
    ::PSY.Line,
) where {T <: ACCEPTED_REAL_TYPES}
    mdl_branch_ode!(
        device_states,
        output_ode,
        voltage_r_from,
        voltage_i_from,
        voltage_r_to,
        voltage_i_to,
        current_r_from,
        current_i_from,
        current_r_to,
        current_i_to,
        branch,
    )

    return
end

function branch!(
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r_from::T,
    voltage_i_from::T,
    voltage_r_to::T,
    voltage_i_to::T,
    current_r_from::AbstractArray{T},
    current_i_from::AbstractArray{T},
    current_r_to::AbstractArray{T},
    current_i_to::AbstractArray{T},
    branch::BranchWrapper,
    ::PSY.Transformer2W,
) where {T <: ACCEPTED_REAL_TYPES}
    mdl_transformer_Lshape_ode!(
        device_states,
        output_ode,
        voltage_r_from,
        voltage_i_from,
        voltage_r_to,
        voltage_i_to,
        current_r_from,
        current_i_from,
        current_r_to,
        current_i_to,
        branch,
    )

    return
end
