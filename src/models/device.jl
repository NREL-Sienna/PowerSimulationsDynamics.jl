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
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    inner_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{DynG},
    t::Float64,
) where {DynG <: PSY.DynamicGenerator, T <: Real}
    inner_vars[VR_gen_var] = voltage_r
    inner_vars[VI_gen_var] = voltage_i

    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_pss_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for AVR
    mdl_avr_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain ODEs for Machine
    mdl_machine_ode!(
        device_states,
        output_ode,
        inner_vars,
        current_r,
        current_i,
        dynamic_device,
    )

    #Obtain ODEs for PSY.Shaft
    mdl_shaft_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)
    return
end

function device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::StaticWrapper{PSY.Source, U},
    t::Float64,
) where {T <: Real, U <: BusCategory}
    mdl_source!(voltage_r, voltage_i, current_r, current_i, device)
    return
end

function device!(
    voltage_r::T,
    voltage_i::T,
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    ::AbstractArray{T},
    device::StaticWrapper{PSY.PowerLoad, U},
    t,
) where {T <: Real, U <: LoadCategory}
    # For now model all loads as constant impedance
    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device.device)
    return
end

function device_mass_matrix_entries!(
    mass_matrix::AbstractArray{Float64},
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
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::AbstractArray{T},
    voltage_i::AbstractArray{T},
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    global_vars::AbstractArray{T},
    dynamic_device::DynamicWrapper{DynI},
    inner_vars::AbstractArray{T},
    t::Float64,
) where {DynI <: PSY.DynamicInverter, T <: Real}
    bus_ix = get_bus_ix(dynamic_device)

    #Obtain global vars
    sys_ω = global_vars[GLOBAL_VAR_SYS_FREQ_INDEX]

    #Update Voltage data
    inner_vars[Vr_inv_var] = voltage_r[bus_ix]
    inner_vars[Vi_inv_var] = voltage_i[bus_ix]

    #Update V_ref
    get_V_ref(dynamic_device)
    inner_vars[V_oc_var] = V_ref

    #Obtain ODES for DC side
    mdl_DCside_ode!(device_states, output_ode, sys_ω, inner_vars, dynamic_device)

    #Obtain ODEs for PLL
    mdl_freq_estimator_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain ODEs for OuterLoop
    mdl_outer_ode!(device_states, output_ode, inner_vars, sys_ω, dynamic_device)

    #Obtain inner controller ODEs and modulation commands
    mdl_inner_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain converter relations
    mdl_converter_ode!(device_states, output_ode, inner_vars, dynamic_device)

    #Obtain ODEs for output filter
    mdl_filter_ode!(
        device_states,
        output_ode,
        current_r,
        current_i,
        inner_vars,
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
    device_states::AbstractArray{T},
    output_ode::AbstractArray{T},
    voltage_r::AbstractArray{T},
    voltage_i::AbstractArray{T},
    current_r::AbstractArray{T},
    current_i::AbstractArray{T},
    ::AbstractArray{T},
    dynamic_device::PSY.PeriodicVariableSource,
    ::AbstractArray{T},
    t::Float64,
) where {T <: Real}
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
    V_R = device_states[1] * cos(device_states[2])
    V_I = device_states[1] * sin(device_states[2])
    output_ode[1] = dV
    output_ode[2] = dθ

    #update current
    R_th = PSY.get_R_th(dynamic_device)
    X_th = PSY.get_X_th(dynamic_device)
    Zmag = R_th^2 + X_th^2
    current_r[1] += R_th * (V_R - voltage_r[1]) / Zmag + X_th * (V_I - voltage_i[1]) / Zmag #in system pu flowing out
    current_i[1] += R_th * (V_I - voltage_i[1]) / Zmag - X_th * (V_R - voltage_r[1]) / Zmag #in system pu flowing out

    return
end
