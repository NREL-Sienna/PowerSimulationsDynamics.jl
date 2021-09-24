function mass_matrix_filter_entries!(
    mass_matrix,
    filt::F,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {F <: PSY.Filter}
    @debug "Using default mass matrix entries $F"
end

function mdl_filter_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter}},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}

    #external_ix = get_input_port_ix(dynamic_device, PSY.LCLFilter)
    #θ_oc = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[Vr_inv_var]
    V_tI = inner_vars[Vi_inv_var]
    Vr_cnv = inner_vars[Vr_cnv_var]
    Vi_cnv = inner_vars[Vi_cnv_var]

    #Get parameters
    filter = PSY.get_filter(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0
    lf = PSY.get_lf(filter)
    rf = PSY.get_rf(filter)
    cf = PSY.get_cf(filter)
    lg = PSY.get_lg(filter)
    rg = PSY.get_rg(filter)
    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = get_system_base_power(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.LCLFilter)

    #Define internal states for filter
    internal_states = @view device_states[local_ix]
    Ir_cnv = internal_states[1]
    Ii_cnv = internal_states[2]
    Vr_filter = internal_states[3]
    Vi_filter = internal_states[4]
    Ir_filter = internal_states[5]
    Ii_filter = internal_states[6]

    #Inputs (control signals) - N/A

    #Compute 6 states ODEs (D'Arco EPSR122 Model)
    #Inverter Output Inductor (internal state)
    #𝜕id_c/𝜕t
    output_ode[local_ix[1]] = (
        ωb / lf * Vr_cnv - ωb / lf * Vr_filter - ωb * rf / lf * Ir_cnv +
        ωb * ω_sys * Ii_cnv
    )
    #𝜕iq_c/𝜕t
    output_ode[local_ix[2]] = (
        ωb / lf * Vi_cnv - ωb / lf * Vi_filter - ωb * rf / lf * Ii_cnv -
        ωb * ω_sys * Ir_cnv
    )
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[local_ix[3]] =
        (ωb / cf * Ir_cnv - ωb / cf * Ir_filter + ωb * ω_sys * Vi_filter)
    #𝜕vq_o/𝜕t
    output_ode[local_ix[4]] =
        (ωb / cf * Ii_cnv - ωb / cf * Ii_filter - ωb * ω_sys * Vr_filter)
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[local_ix[5]] = (
        ωb / lg * Vr_filter - ωb / lg * V_tR - ωb * rg / lg * Ir_filter +
        ωb * ω_sys * Ii_filter
    )
    #𝜕iq_o/𝜕t
    output_ode[local_ix[6]] = (
        ωb / lg * Vi_filter - ωb / lg * V_tI - ωb * rg / lg * Ii_filter -
        ωb * ω_sys * Ir_filter
    )

    #Update inner_vars
    inner_vars[Vr_filter_var] = Vr_filter
    inner_vars[Vi_filter_var] = Vi_filter

    #Compute current from the inverter to the grid
    I_RI = (basepower / sys_Sbase) * [Ir_filter; Ii_filter]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]
end

function mdl_filter_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    inner_vars,
    ω_sys,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.RLFilter}},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
}
    #Obtain inner variables for component
    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = get_system_base_power(dynamic_device)
    ratio_power = basepower / sys_Sbase

    #Obtain parameters
    filt = PSY.get_filter(dynamic_device)
    rf = PSY.get_rf(filt)
    lf = PSY.get_lf(filt)

    #Compute output currents
    if lf != 0.0 || rf != 0.0
        Vr_cnv = inner_vars[Vr_cnv_var]
        Vi_cnv = inner_vars[Vi_cnv_var]
        Vr_inv = inner_vars[Vr_inv_var]
        Vi_inv = inner_vars[Vi_inv_var]
        Zmag_squared = rf^2 + lf^2
        Ir_filt = (1.0 / Zmag_squared) * ((Vr_cnv - Vr_inv) * rf + (Vi_cnv - Vi_inv) * lf)
        Ii_filt = (1.0 / Zmag_squared) * ((Vi_cnv - Vi_inv) * rf - (Vr_cnv - Vr_inv) * lf)
    else
        Ir_filt = inner_vars[Ir_cnv_var]
        Ii_filt = inner_vars[Ii_cnv_var]
    end

    #Update Inner Vars
    inner_vars[Ir_inv_var] = Ir_filt
    inner_vars[Ii_inv_var] = Ii_filt

    #Update current
    current_r[1] += ratio_power * Ir_filt
    current_i[1] += ratio_power * Ii_filt
end
