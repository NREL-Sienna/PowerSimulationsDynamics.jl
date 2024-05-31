function mass_matrix_DCside_entries!(
    mass_matrix,
    dc_side::DC,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {DC <: PSY.DCSource}
    CRC.@ignore_derivatives @debug "Using default mass matrix entries $DC"
end

function mdl_DCside_ode!(
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ::AbstractArray{<:ACCEPTED_REAL_TYPES},
    p::AbstractArray{<:ACCEPTED_REAL_TYPES},
    ω_sys::ACCEPTED_REAL_TYPES,
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicInverter{C, O, IC, PSY.FixedDCSource, P, F, L},
    },
    h,
    t,
) where {
    C <: PSY.Converter,
    O <: PSY.OuterControl,
    IC <: PSY.InnerControl,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
    L <: Union{Nothing, PSY.InverterLimiter},
}
    #Update inner_vars
    inner_vars[Vdc_var] = p[:params][:DCSource][:voltage]
end
