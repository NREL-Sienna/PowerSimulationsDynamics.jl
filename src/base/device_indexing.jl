function _attach_control_refs!(device::PSY.StaticInjection)
    dynamic_device = PSY.get_dynamic_injector(device)
    dynamic_device.ext[CONTROL_REFS] = [
        PSY.get_V_ref(dynamic_device),
        PSY.get_ω_ref(dynamic_device),
        PSY.get_P_ref(dynamic_device),
        PSY.get_reactive_power(device),
    ]
    return
end

function _attach_control_refs!(device::PSY.Source)
    device.ext[CONTROL_REFS] =
        [PSY.get_internal_voltage(device), PSY.get_internal_angle(device)]
    return
end

function _attach_control_refs!(device::PSY.ElectricLoad)
    return
end

function index_static_injection(inputs::SimulationInputs, ::PSY.StaticInjection)
    return
end
