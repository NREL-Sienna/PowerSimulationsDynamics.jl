function initialize_avr!(device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = get_inner_vars(device)[Vf_var]
    PSY.set_V_ref!(PSY.get_avr(device), Vf)
    #Update Control Refs
    device.ext[CONTROL_REFS][1] = Vf
end