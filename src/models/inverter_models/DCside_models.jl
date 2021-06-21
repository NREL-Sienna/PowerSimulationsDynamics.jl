function mass_matrix_DCside_entries!(
    mass_matrix,
    dc_side::DC,
    global_index::Dict{Symbol, Int64},
) where {DC <: PSY.DCSource}
    @debug "Using default mass matrix entries $DC"
end

function mdl_DCside_ode!(
    device_states,
    output_ode,
    f0,
    ω_sys,
    dynamic_device::PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F},
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Update inner_vars
    set_inner_vars!(
        dynamic_device,
        Vdc_var,
        PSY.get_voltage(PSY.get_dc_source(dynamic_device)),
    )
end
