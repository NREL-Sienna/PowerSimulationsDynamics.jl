function mass_matrix_filter_entries!(
    mass_matrix,
    filt::F,
    global_index::Base.ImmutableDict{Symbol, Int64},
    f0::Float64,
) where {F <: PSY.Filter}
    @debug "Using default mass matrix entries $F"
    return
end

function mass_matrix_filter_entries!(
    mass_matrix,
    filt::PSY.LCLFilter,
    global_index::Base.ImmutableDict{Symbol, Int64},
    f0::Float64,
)
    ext = PSY.get_ext(filt)
    bool_val = get(ext, "is_filter_differential", 1.0)
    ωb = 2 * pi * f0
    mass_matrix[global_index[:ir_cnv], global_index[:ir_cnv]] =
        bool_val * PSY.get_lf(filt) / ωb
    mass_matrix[global_index[:ii_cnv], global_index[:ii_cnv]] =
        bool_val * PSY.get_lf(filt) / ωb
    mass_matrix[global_index[:vr_filter], global_index[:vr_filter]] =
        bool_val * PSY.get_cf(filt) / ωb
    mass_matrix[global_index[:vi_filter], global_index[:vi_filter]] =
        bool_val * PSY.get_cf(filt) / ωb
    mass_matrix[global_index[:ir_filter], global_index[:ir_filter]] =
        bool_val * PSY.get_lg(filt) / ωb
    mass_matrix[global_index[:ii_filter], global_index[:ii_filter]] =
        bool_val * PSY.get_lg(filt) / ωb
    return
end

function mdl_filter_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.LCLFilter, L}},
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}

    #external_ix = get_input_port_ix(dynamic_device, PSY.LCLFilter)
    #θ_oc = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[Vr_inv_var]
    V_tI = inner_vars[Vi_inv_var]
    Vr_cnv = inner_vars[Vr_cnv_var]
    Vi_cnv = inner_vars[Vi_cnv_var]

    #Get parameters
    f0 = get_system_base_frequency(dynamic_device)
    ωb = 2 * pi * f0
    lf = p[:params][:Filter][:lf]
    rf = p[:params][:Filter][:rf]
    cf = p[:params][:Filter][:cf]
    lg = p[:params][:Filter][:lg]
    rg = p[:params][:Filter][:rg]
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
    output_ode[local_ix[1]] = (Vr_cnv - Vr_filter - rf * Ir_cnv + lf * ω_sys * Ii_cnv)
    #𝜕iq_c/𝜕t
    output_ode[local_ix[2]] = (Vi_cnv - Vi_filter - rf * Ii_cnv - lf * ω_sys * Ir_cnv)
    #LCL Capacitor (internal state)
    #𝜕vd_o/𝜕t
    output_ode[local_ix[3]] = (Ir_cnv - Ir_filter + cf * ω_sys * Vi_filter)
    #𝜕vq_o/𝜕t
    output_ode[local_ix[4]] = (Ii_cnv - Ii_filter - cf * ω_sys * Vr_filter)
    #Grid Inductance (internal state)
    #𝜕id_o/𝜕t
    output_ode[local_ix[5]] = (Vr_filter - V_tR - rg * Ir_filter + lg * ω_sys * Ii_filter)
    #𝜕iq_o/𝜕t
    output_ode[local_ix[6]] = (Vi_filter - V_tI - rg * Ii_filter - lg * ω_sys * Ir_filter)

    #Update inner_vars
    inner_vars[Vr_filter_var] = Vr_filter
    inner_vars[Vi_filter_var] = Vi_filter

    #Compute current from the inverter to the grid
    I_RI = (basepower / sys_Sbase) * [Ir_filter; Ii_filter]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]
    return
end

function mdl_filter_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, DC, P, PSY.RLFilter, L}},
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    L <: Union{Nothing, PSY.OutputCurrentLimiter},
}
    #Obtain inner variables for component
    basepower = PSY.get_base_power(dynamic_device)
    sys_Sbase = get_system_base_power(dynamic_device)
    ratio_power = basepower / sys_Sbase

    #Obtain parameters
    rf = p[:params][:Filter][:rf]
    lf = p[:params][:Filter][:lf]

    Vr_cnv = inner_vars[Vr_cnv_var]
    Vi_cnv = inner_vars[Vi_cnv_var]
    Vr_inv = inner_vars[Vr_inv_var]
    Vi_inv = inner_vars[Vi_inv_var]

    #Compute output currents
    if lf != 0.0 || rf != 0.0
        Zmag_squared = rf^2 + lf^2
        Ir_filt = (1.0 / Zmag_squared) * ((Vr_cnv - Vr_inv) * rf + (Vi_cnv - Vi_inv) * lf)
        Ii_filt = (1.0 / Zmag_squared) * ((Vi_cnv - Vi_inv) * rf - (Vr_cnv - Vr_inv) * lf)
    else
        #Obtain converter
        IS.@assert_op Vr_cnv == Vr_inv
        IS.@assert_op Vi_cnv == Vi_inv
        converter = PSY.get_converter(dynamic_device)
        R_source = PSY.get_R_source(converter)
        X_source = PSY.get_X_source(converter)
        Z_source_sq = R_source^2 + X_source^2
        I_aux_r = (R_source * Vr_cnv + X_source * Vi_cnv) / Z_source_sq
        I_aux_i = (R_source * Vi_cnv - X_source * Vr_cnv) / Z_source_sq
        Ir_filt = inner_vars[Ir_cnv_var] - I_aux_r
        Ii_filt = inner_vars[Ii_cnv_var] - I_aux_i
    end
    #Update Inner Vars
    inner_vars[Ir_inv_var] = Ir_filt
    inner_vars[Ii_inv_var] = Ii_filt

    #Update current
    current_r[1] += ratio_power * Ir_filt
    current_i[1] += ratio_power * Ii_filt
    return
end
