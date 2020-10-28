function device!(
    x,
    output_ode,
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int64},
    ode_range::UnitRange{Int64},
    dyn_data::DynG,
    inputs::SimulationInputs,
) where {DynG <: PSY.DynamicGenerator}
    sys = get_system(inputs)
    #Obtain local device states
    n_states = PSY.get_n_states(dyn_data)
    device_states = @view x[ix_range]

    #Obtain references
    sys_Sbase = PSY.get_base_power(sys)
    sys_f0 = PSY.get_frequency(sys)
    sys_ω = get_ω_sys(inputs)

    #Update Voltage data
    get_inner_vars(dyn_data)[VR_gen_var] = voltage_r[1]
    get_inner_vars(dyn_data)[VI_gen_var] = voltage_i[1]

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, view(output_ode, ode_range), dyn_data)

    #Obtain ODEs for AVR
    mdl_pss_ode!(device_states, view(output_ode, ode_range), sys_ω, dyn_data)

    #Obtain ODEs for AVR
    mdl_avr_ode!(device_states, view(output_ode, ode_range), dyn_data)

    #Obtain ODEs for Machine
    mdl_machine_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        sys_Sbase,
        sys_f0,
        dyn_data,
    )

    #Obtain ODEs for PSY.Shaft
    mdl_shaft_ode!(device_states, view(output_ode, ode_range), sys_f0, sys_ω, dyn_data)

    return
end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.Source,
    inputs::SimulationInputs,
)
    sys = get_system(inputs)
    mdl_source!(voltage_r, voltage_i, current_r, current_i, device, sys)
    return
end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.PowerLoad,
    inputs::SimulationInputs,
)
    sys = get_system(inputs)
    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device, sys)
    return
end

function device!(
    x,
    output_ode::Vector{T},
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int64},
    ode_range::UnitRange{Int64},
    dyn_data::DynI,
    inputs::SimulationInputs,
) where {DynI <: PSY.DynamicInverter, T <: Real}
    sys = get_system(inputs)
    #Obtain local device states
    n_states = PSY.get_n_states(dyn_data)
    device_states = @view x[ix_range]

    #Obtain references
    Sbase = PSY.get_base_power(sys)
    sys_f0 = PSY.get_frequency(sys)
    sys_ω = get_ω_sys(inputs)

    #Update Voltage data
    get_inner_vars(dyn_data)[VR_inv_var] = voltage_r[1]
    get_inner_vars(dyn_data)[VI_inv_var] = voltage_i[1]

    #Update V_ref
    V_ref = PSY.get_ext(dyn_data)[CONTROL_REFS][V_ref_index]
    get_inner_vars(dyn_data)[V_oc_var] = V_ref

    #Obtain ODES for DC side
    mdl_DCside_ode!(dyn_data)

    #Obtain ODEs for PLL
    mdl_freq_estimator_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        dyn_data,
    )

    #Obtain ODEs for OuterLoop
    mdl_outer_ode!(device_states, view(output_ode, ode_range), sys_f0, sys_ω, dyn_data)

    #Obtain inner controller ODEs and modulation commands
    mdl_inner_ode!(device_states, view(output_ode, ode_range), dyn_data)

    #Obtain converter relations
    mdl_converter_ode!(dyn_data)

    #Obtain ODEs for output filter
    mdl_filter_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        Sbase,
        sys_f0,
        sys_ω,
        dyn_data,
    )

    return
end
