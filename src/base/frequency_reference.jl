struct FixedFrequency end

struct ReferenceBus end

function get_frequency_reference!(
    ::Type{FixedFrequency},
    wrapped_injectors,
    static_injection_data,
)
    ref_devices = filter(
        x -> PSY.get_bustype(PSY.get_bus(x)) == PSY.BusTypes.REF,
        static_injection_data,
    )

    if length(ref_devices) < 1
        throw(
            IS.ConflictingInputsError(
                "InfiniteBus model requires at least one bus of type BusTypes.REF with a Source connected to it",
            ),
        )
    end

    return 0
end

function get_frequency_reference!(::Type{ReferenceBus}, wrapped_injectors)
    reference = -1
    ref_devices = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x ->
            PSY.get_bustype(PSY.get_bus(x)) == PSY.BusTypes.REF &&
                !isa(x, PSY.ElectricLoad),
    )
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
        reference = get_frequency_reference!(inputs, first(ref_devices))
    end
    return
end
