function mass_matrix_DCside_entries!(
    mass_matrix,
    dc_side::DC,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {DC <: PSY.DCSource}
    @debug "Using default mass matrix entries $DC"
end

function mdl_DCside_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F}},
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Update inner_vars
    inner_vars[Vdc_var] = PSY.get_voltage(PSY.get_dc_source(dynamic_device))
end
