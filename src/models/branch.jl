function branch!(
    x,
    dx,
    output_ode::Vector{T},
    V_r_from,
    V_i_from,
    V_r_to,
    V_i_to,
    current_r_from,
    current_i_from,
    current_r_to,
    current_i_to,
    ix_range::UnitRange{Int},
    ode_range::UnitRange{Int},
    branch::PSY.DynamicBranch,
    inputs::SimulationInputs,
) where {T <: Real}
    #Obtain local device states
    device_states = @view x[ix_range]

    #Obtain references
    sys_f = get_base_frequency(inputs)

    mdl_line_ode!(
        device_states,
        view(output_ode, ode_range),
        V_r_from,
        V_i_from,
        V_r_to,
        V_i_to,
        current_r_from,
        current_i_from,
        current_r_to,
        current_i_to,
        sys_f,
        branch,
    )

    return
end
