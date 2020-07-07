function initialize_tg!(device_states,
    device::PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    tg = PSY.get_prime_mover(device)
    τm0 = get_inner_vars(device)[τm_var]
    eff = PSY.get_efficiency(tg)
    P_ref = τm0 / eff
    PSY.set_P_ref!(tg, P_ref)
    #Update Control Refs
    device.ext[CONTROL_REFS][P_ref_index] = P_ref
end