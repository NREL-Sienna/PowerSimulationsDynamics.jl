function device_mass_matrix_entries!(::AbstractArray, ::DynamicWrapper{T}) where {T}
    error("Mass Matrix not implemented for models $T")
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{DynG},
) where {DynG <: PSY.DynamicGenerator}
    global_index = get_global_index(dynamic_device)
    mass_matrix_tg_entries!(mass_matrix, PSY.get_prime_mover(dynamic_device), global_index)
    mass_matrix_pss_entries!(mass_matrix, PSY.get_pss(dynamic_device), global_index)
    mass_matrix_avr_entries!(mass_matrix, PSY.get_avr(dynamic_device), global_index)
    mass_matrix_machine_entries!(mass_matrix, PSY.get_machine(dynamic_device), global_index)
    mass_matrix_shaft_entries!(mass_matrix, PSY.get_shaft(dynamic_device), global_index)
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
    dynamic_device::DynamicWrapper{DynG},
    inputs::SimulationInputs,
    t,
) where {DynG <: PSY.DynamicGenerator, T <: Real}
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
    t,
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
    t,
)
    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray,
    dynamic_device::DynamicWrapper{DynI},
) where {DynI <: PSY.DynamicInverter}
    global_index = get_global_index(dynamic_device)
    mass_matrix_DCside_entries!(
        mass_matrix,
        PSY.get_dc_source(dynamic_device),
        global_index,
    )
    mass_matrix_freq_estimator_entries!(
        mass_matrix,
        PSY.get_freq_estimator(dynamic_device),
        global_index,
    )
    mass_matrix_outer_entries!(
        mass_matrix,
        PSY.get_outer_control(dynamic_device),
        global_index,
    )
    mass_matrix_inner_entries!(
        mass_matrix,
        PSY.get_inner_control(dynamic_device),
        global_index,
    )
    mass_matrix_converter_entries!(
        mass_matrix,
        PSY.get_converter(dynamic_device),
        global_index,
    )
    mass_matrix_filter_entries!(mass_matrix, PSY.get_filter(dynamic_device), global_index)
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
    dynamic_device::DynamicWrapper{DynI},
    inputs::SimulationInputs,
    cache::Cache,
    t,
) where {DynI <: PSY.DynamicInverter, T <: Real}
    #Obtain local device states
    device_states = @view x[ix_range]

    #Obtain references
    Sbase = get_base_power(inputs)
    sys_f0 = get_base_frequency(inputs)
    sys_ω = get_ω_sys(inputs)

    #Update Voltage data
    get_inner_vars(dynamic_device)[Vr_inv_var] = voltage_r[1]
    get_inner_vars(dynamic_device)[Vi_inv_var] = voltage_i[1]

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
    mdl_converter_ode!(device_states, view(output_ode, ode_range), dynamic_device)

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

function device_mass_matrix_entries!(
    sim_inputs::SimulationInputs,
    dynamic_device::PSY.PeriodicVariableSource,
)
    mass_matrix = get_mass_matrix(sim_inputs)
    global_index = get_global_index(sim_inputs)[PSY.get_name(dynamic_device)]
    mass_matrix_pvs_entries!(mass_matrix, dynamic_device, global_index)
    return
end

function mass_matrix_pvs_entries!(
    mass_matrix,
    pvs::PSY.PeriodicVariableSource,
    global_index::Dict{Symbol, Int64},
)
    @debug "Using default mass matrix entries $pvs"
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
    dynamic_device::PSY.PeriodicVariableSource,
    inputs::SimulationInputs,
    t,
) where {T <: Real}
    #Obtain local device states
    internal_states = @view x[ix_range]
    internal_ode = @view output_ode[ode_range]
    ω_θ = PSY.get_internal_angle_frequencies(dynamic_device)
    ω_V = PSY.get_internal_angle_frequencies(dynamic_device)

    dV = 0
    for (ix, A) in enumerate(PSY.get_internal_voltage_coefficients(dynamic_device))
        t <= 0 && continue
        dV += ω_V[ix] * (A[1] * cos(ω_V[ix] * t) - A[2] * sin(ω_V[ix] * t))
    end

    dθ = 0
    for (ix, A) in enumerate(PSY.get_internal_angle_coefficients(dynamic_device))
        t <= 0 && continue
        dθ += ω_θ[ix] * (A[1] * cos(ω_θ[ix] * t) - A[2] * sin(ω_θ[ix] * t))
    end

    # Internal Voltage states
    V_R = internal_states[1] * cos(internal_states[2])
    V_I = internal_states[1] * sin(internal_states[2])
    internal_ode[1] = dV
    internal_ode[2] = dθ

    #update current
    R_th = PSY.get_R_th(dynamic_device)
    X_th = PSY.get_X_th(dynamic_device)
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end
