function device!(
    x,
    output_ode::Vector{T},
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    ix_range::UnitRange{Int64},
    ode_range::UnitRange{Int64},
    device::DynG,
    sys::PSY.System,
) where {DynG <: PSY.DynamicGenerator, T <: Real}
    #Obtain local device states
    n_states = PSY.get_n_states(device)
    device_states = @view x[ix_range]

    #Obtain references
    sys_Sbase = PSY.get_basepower(sys)
    sys_f0 = PSY.get_frequency(sys)
    sys_ω = get_ω_sys(sys)

    #Update Voltage data
    get_inner_vars(device)[VR_gen_var] = voltage_r[1]
    get_inner_vars(device)[VI_gen_var] = voltage_i[1]

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, view(output_ode, ode_range), device)

    #Obtain ODEs for AVR
    mdl_pss_ode!(device_states, view(output_ode, ode_range), sys_ω, device)

    #Obtain ODEs for AVR
    mdl_avr_ode!(device_states, view(output_ode, ode_range), device)

    #Obtain ODEs for Machine
    mdl_machine_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        sys_Sbase,
        sys_f0,
        device,
    )

    #Obtain ODEs for PSY.Shaft
    mdl_shaft_ode!(device_states, view(output_ode, ode_range), sys_f0, sys_ω, device)

    return

end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.Source,
    sys::PSY.System,
)

    mdl_source!(voltage_r, voltage_i, current_r, current_i, device, sys)

    return
end

function device!(
    voltage_r,
    voltage_i,
    current_r,
    current_i,
    device::PSY.PowerLoad,
    sys::PSY.System,
)

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
    device::DynI,
    sys::PSY.System,
) where {DynI <: PSY.DynamicInverter, T <: Real}

    #Obtain local device states
    n_states = PSY.get_n_states(device)
    device_states = @view x[ix_range]

    #Obtain references
    Sbase = PSY.get_basepower(sys)
    sys_f0 = PSY.get_frequency(sys)
    sys_ω = get_ω_sys(sys)

    #Update Voltage data
    get_inner_vars(device)[VR_inv_var] = voltage_r[1]
    get_inner_vars(device)[VI_inv_var] = voltage_i[1]

    #Update V_ref
    V_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]
    get_inner_vars(device)[V_oc_var] = V_ref

    #Obtain ODES for DC side
    mdl_DCside_ode!(device)

    #Obtain ODEs for PLL
    mdl_freq_estimator_ode!(
        device_states,
        view(output_ode, ode_range),
        sys_f0,
        sys_ω,
        device,
    )

    #Obtain ODEs for OuterLoop
    mdl_outer_ode!(device_states, view(output_ode, ode_range), sys_f0, sys_ω, device)

    #Obtain inner controller ODEs and modulation commands
    mdl_inner_ode!(device_states, view(output_ode, ode_range), device)

    #Obtain converter relations
    mdl_converter_ode!(device)

    #Obtain ODEs for output filter
    mdl_filter_ode!(
        device_states,
        view(output_ode, ode_range),
        current_r,
        current_i,
        Sbase,
        sys_f0,
        sys_ω,
        device,
    )

    return
end
