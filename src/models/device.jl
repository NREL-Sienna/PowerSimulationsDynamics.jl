function device!(
    x,
    output_ode,
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int},
    ode_range::UnitRange{Int},
    dynamic_device::DynG,
    inputs::SimulationInputs,
) where {DynG <: PSY.DynamicGenerator}
    #Obtain local device states
    device_states = @view x[ix_range]

    #Obtain references
    sys_Sbase = get_base_power(inputs)
    sys_f0 = get_base_frequency(inputs)
    sys_ω = get_ω_sys(inputs)

    #Update Voltage data
    get_inner_vars(dynamic_device)[VR_gen_var] = voltage_r[1]
    get_inner_vars(dynamic_device)[VI_gen_var] = voltage_i[1]

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, view(output_ode, ode_range), sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_pss_ode!(device_states, view(output_ode, ode_range), sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_avr_ode!(device_states, view(output_ode, ode_range), dynamic_device)

    #Obtain ODEs for Machine
    mdl_machine_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        sys_Sbase,
        sys_f0,
        dynamic_device,
    )

    #Obtain ODEs for PSY.Shaft
    mdl_shaft_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        dynamic_device,
    )

    return
end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.Source,
    ::SimulationInputs,
)
    mdl_source!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.PowerLoad,
    ::SimulationInputs,
)
    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function device!(
    x,
    output_ode::Vector{T},
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int},
    ode_range::UnitRange{Int},
    dynamic_device::DynI,
    inputs::SimulationInputs,
) where {DynI <: PSY.DynamicInverter, T <: Real}
    #Obtain local device states
    device_states = @view x[ix_range]

    #Obtain references
    Sbase = get_base_power(inputs)
    sys_f0 = get_base_frequency(inputs)
    sys_ω = get_ω_sys(inputs)

    #Update Voltage data
    get_inner_vars(dynamic_device)[VR_inv_var] = voltage_r[1]
    get_inner_vars(dynamic_device)[VI_inv_var] = voltage_i[1]

    #Update V_ref
    V_ref = PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index]
    get_inner_vars(dynamic_device)[V_oc_var] = V_ref

    #Obtain ODES for DC side
    mdl_DCside_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        dynamic_device,
    )

    #Obtain ODEs for PLL
    mdl_freq_estimator_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        dynamic_device,
    )

    #Obtain ODEs for OuterLoop
    mdl_outer_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        dynamic_device,
    )

    #Obtain inner controller ODEs and modulation commands
    mdl_inner_ode!(device_states, view(output_ode, ode_range), dynamic_device)

    #Obtain converter relations
    mdl_converter_ode!(dynamic_device)

    #Obtain ODEs for output filter
    mdl_filter_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        Sbase,
        sys_f0,
        sys_ω,
        dynamic_device,
    )

    return
end
