function initialize_converter!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{PSY.AverageConverter, O, IC, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.InverterLimiter},
} end

function initialize_converter!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{PSY.RenewableEnergyConverterTypeA, O, IC, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.InverterLimiter},
}
    #Get inner vars
    V_R = inner_vars[Vr_cnv_var]
    V_I = inner_vars[Vi_cnv_var]
    θ = atan(V_I, V_R)
    bus = PSY.get_bus(static)
    bus_angle = PSY.get_angle(bus)
    @assert isapprox(θ, bus_angle)
    V_t = sqrt(V_R^2 + V_I^2)
    Iq_external = inner_vars[Ii_cnv_var]
    Ip_external = inner_vars[Ir_cnv_var]
    #Reference Transformation
    Ip = Ip_external * cos(-θ) - Iq_external * sin(-θ)
    Iq = Ip_external * sin(-θ) + Iq_external * cos(-θ)

    #Get Converter parameters
    converter = PSY.get_converter(dynamic_device)
    local_ix_params =
        get_local_parameter_ix(dynamic_device, PSY.RenewableEnergyConverterTypeA)
    internal_params = @view device_parameters[local_ix_params]
    T_g,
    Rrpwr,
    Brkpt,
    Zerox,
    Lvpl1,
    Vo_lim,
    Lv_pnt0,
    Lv_pnt1,
    Io_lim,
    T_fltr,
    K_hv,
    Iqr_min,
    Iqr_max,
    Accel,
    Q_ref,
    R_source,
    X_source = internal_params

    # Lv_pnt0 is unused in the initialization
    _, Lv_pnt1 = PSY.get_Lv_pnts(converter)

    if (Iq < Io_lim) || (V_t > Vo_lim) || (V_t < Lv_pnt1)
        CRC.@ignore_derivatives @error(
            "Power flow solution outside of inverter limits $(PSY.get_name(static)). Update parameters."
        )
    end

    #Update converter states
    converter_ix = get_local_state_ix(dynamic_device, PSY.RenewableEnergyConverterTypeA)
    converter_states = @view device_states[converter_ix]
    converter_states[1] = Ip #Ip_cnv
    converter_states[2] = -Iq #Iq_cnv
    converter_states[3] = V_t #Vmeas

    #Update entry to converter currents (inner currents)
    inner_vars[Id_ic_var] = Ip
    # Multiplied by -1 because PSS/e doesn't do a conjugate.
    inner_vars[Iq_ic_var] = -Iq
    return
end

function initialize_converter!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{PSY.RenewableEnergyVoltageConverterTypeA, O, IC, DC, P, F, L},
    },
    inner_vars::AbstractVector,
) where {
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.InverterLimiter},
}
    #Get inner vars
    Vr_filter = inner_vars[Vr_filter_var]
    Vi_filter = inner_vars[Vi_filter_var]
    θ = atan(Vi_filter, Vr_filter)
    V_t = sqrt(Vr_filter^2 + Vi_filter^2)
    Iq_external = inner_vars[Ii_cnv_var]
    Ip_external = inner_vars[Ir_cnv_var]
    #Reference Transformation
    Ip = Ip_external * cos(-θ) - Iq_external * sin(-θ)
    Iq = Ip_external * sin(-θ) + Iq_external * cos(-θ)
    converter = PSY.get_converter(dynamic_device)
    Io_lim = PSY.get_Io_lim(converter)
    Vo_lim = PSY.get_Vo_lim(converter)

    # Lv_pnt0 is unused in the initialization
    _, Lv_pnt1 = PSY.get_Lv_pnts(converter)

    if (Iq < Io_lim) || (V_t > Vo_lim) || (V_t < Lv_pnt1)
        CRC.@ignore_derivatives @error(
            "Power flow solution outside of inverter limits $(PSY.get_name(static)). Update parameters."
        )
    end

    #Update converter states
    converter_ix =
        get_local_state_ix(dynamic_device, PSY.RenewableEnergyVoltageConverterTypeA)
    converter_states = @view device_states[converter_ix]
    converter_states[1] = Ip #Ip_cnv
    converter_states[2] = -Iq #Iq_cnv
    converter_states[3] = V_t #Vmeas

    #Update entry to converter currents (inner currents)
    inner_vars[Id_ic_var] = Ip
    # Multiplied by -1 because PSS/e doesn't do a conjugate.
    inner_vars[Iq_ic_var] = -Iq
    return
end
