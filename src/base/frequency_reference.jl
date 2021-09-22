struct FixedFrequency end
struct ReferenceBus end

function _get_frequency_state(d::DynamicWrapper{T}) where {T <: PSY.DynamicGenerator}
    return get_global_index(d)[:ω]
end

function _get_frequency_state(d::DynamicWrapper{T}) where {T <: PSY.DynamicInverter}
    return get_global_index(d)[:ω_oc]
end

function _get_frequency_state(d::DynamicWrapper{PSY.PeriodicVariableSource})
    return 0
end

function get_frequency_reference!(
    ::Type{FixedFrequency},
    ::Vector,
    static_injection_data::Vector,
)
    ref_devices = filter(x -> get_bus_category(x) == SLACKBus, static_injection_data)
    if length(ref_devices) < 1
        throw(
            IS.ConflictingInputsError(
                "InfiniteBus model requires at least one bus of type BusTypes.REF with a Source connected to it",
            ),
        )
    end
    return 0
end

function get_frequency_reference(
    ::Type{ReferenceBus},
    dynamic_injection_data::Vector,
    static_injection_data::Vector,
)
    reference = -1
    static_ref_devices = filter(x -> get_bus_category(x) == SLACKBus, static_injection_data)
    if isempty(static_ref_devices)
        ref_devices = filter(x -> get_bus_category(x) == SLACKBus, dynamic_injection_data)
        if length(ref_devices) > 1
            throw(
                IS.ConflictingInputsError(
                    "More than one source or generator in the REF Bus is not supported in ReferenceBus model",
                ),
            )
        elseif length(ref_devices) < 1
            throw(
                IS.ConflictingInputsError(
                    "ReferenceBus model requires at least one bus of type BusTypes.REF with a DynamicInjection or Source device connected to it",
                ),
            )
        else
            reference = _get_frequency_state(first(ref_devices))
        end
    else
        @warn(
            "The reference Bus has a Source connected to it. The frequency reference model will change to FixedFrequency"
        )
        reference = 0
    end
    return reference
end
