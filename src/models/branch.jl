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
    ix_range::UnitRange{Int64},
    ode_range::UnitRange{Int64},
    branch::PSY.DynamicBranch,
    inputs::SimulationInputs,
) where {T <: Real}
    sys = get_system(inputs)
    #Obtain local device states
    n_states = PSY.get_n_states(branch)
    device_states = @view x[ix_range]

    #Obtain references
    Sbase = PSY.get_base_power(sys)
    sys_f = PSY.get_frequency(sys)

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
