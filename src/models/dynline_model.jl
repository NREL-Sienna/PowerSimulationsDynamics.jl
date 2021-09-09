function mdl_line_ode!(
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
) where {T <: Real}
    L = PSY.get_x(branch)
    R = PSY.get_r(branch)
    ω_b = get_system_base_frequency(branch) * 2 * π

    Il_r = device_states[1]
    Il_i = device_states[2]
    output_ode[1] = (ω_b / L) * ((voltage_r_from - voltage_r_to) - (R * Il_r - L * Il_i))
    output_ode[2] = (ω_b / L) * ((voltage_i_from - voltage_i_to) - (R * Il_i + L * Il_r))

    current_r_from[1] -= Il_r
    current_i_from[1] -= Il_i
    current_r_to[1] += Il_r
    current_i_to[1] += Il_i
end
