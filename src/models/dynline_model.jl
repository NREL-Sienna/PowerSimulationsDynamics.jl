function mdl_branch_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{T},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    voltage_r_from::T,
    voltage_i_from::T,
    voltage_r_to::T,
    voltage_i_to::T,
    current_r_from::AbstractArray{T},
    current_i_from::AbstractArray{T},
    current_r_to::AbstractArray{T},
    current_i_to::AbstractArray{T},
    branch::BranchWrapper,
) where {T <: ACCEPTED_REAL_TYPES}
    R = p[:params][:r]
    L = p[:params][:x]

    ω_b = get_system_base_frequency(branch) * 2 * π

    Il_r = device_states[1]
    Il_i = device_states[2]
    output_ode[1] = (ω_b / L) * ((voltage_r_from - voltage_r_to) - (R * Il_r - L * Il_i))
    output_ode[2] = (ω_b / L) * ((voltage_i_from - voltage_i_to) - (R * Il_i + L * Il_r))

    current_r_from[1] -= Il_r
    current_i_from[1] -= Il_i
    current_r_to[1] += Il_r
    current_i_to[1] += Il_i
    return
end

function mdl_transformer_Lshape_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{T},
    device_parameters::AbstractArray{<:ACCEPTED_REAL_TYPES},
    voltage_r_from::T,
    voltage_i_from::T,
    voltage_r_to::T,
    voltage_i_to::T,
    current_r_from::AbstractArray{T},
    current_i_from::AbstractArray{T},
    current_r_to::AbstractArray{T},
    current_i_to::AbstractArray{T},
    branch::BranchWrapper,
) where {T <: ACCEPTED_REAL_TYPES}
    L = PSY.get_x(branch)
    R = PSY.get_r(branch)
    ω_b = get_system_base_frequency(branch) * 2 * π
    dyn_branch = branch.branch
    br = dyn_branch.branch
    shunt = PSY.get_primary_shunt(br)
    Il_r = device_states[1]
    Il_i = device_states[2]
    Ishunt_r = device_states[3]
    Ishunt_i = device_states[3]
    output_ode[1] = (ω_b / L) * ((voltage_r_from - voltage_r_to) - (R * Il_r - L * Il_i))
    output_ode[2] = (ω_b / L) * ((voltage_i_from - voltage_i_to) - (R * Il_i + L * Il_r))
    output_ode[3] = (ω_b / shunt) * voltage_r_from + L * Ishunt_i
    output_ode[4] = (ω_b / shunt) * voltage_i_from - L * Ishunt_r

    current_r_from[1] -= (Il_r + Ishunt_r)
    current_i_from[1] -= (Il_r + Ishunt_i)
    current_r_to[1] += Il_r
    current_i_to[1] += Il_i
    return
end
