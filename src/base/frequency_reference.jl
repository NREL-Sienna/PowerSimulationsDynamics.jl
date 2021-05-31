struct ReferenceBus end

function attach_frequency_reference!(inputs::SimulationInputs, d::PSY.DynamicGenerator)
    inputs.global_vars[:ω_sys_index] = inputs.global_index[PSY.get_name(d)][:ω]
    return
end

function attach_frequency_reference!(inputs::SimulationInputs, d::PSY.DynamicInverter)
    inputs.global_vars[:ω_sys_index] = inputs.global_index[PSY.get_name(d)][:ω_oc]
    return
end

function attach_frequency_reference!(inputs::SimulationInputs, ::PSY.Source)
    inputs.global_vars[:ω_sys_index] = 0
    return
end

function set_frequency_reference!(::Type{ReferenceBus}, inputs::SimulationInputs, sys::PSY.system)
    ref_devices = PSY.get_components(PSY.StaticInjection, sys, x -> PSY.get_bustype(PSY.get_bus(x)) == PSY.BusTypes.REF)
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
        attach_frequency_reference!(inputs, first(ref_devices))
    end
    return
end
