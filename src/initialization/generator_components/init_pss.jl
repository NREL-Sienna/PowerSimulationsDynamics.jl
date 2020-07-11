function initialize_pss!(
    device_states,
    device::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

end

#Currently not working properly.
#=
function initialize_pss!(
    device_states,
    device::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    ω0 = PSY.get_ext(device)[CONTROL_REFS][ω_ref_index]
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen) / PSY.get_basepower(static_gen)
    τe = get_inner_vars(device)[τe_var]

    #Get parameters
    pss = PSY.get_pss(device)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #It assumes that the system is balanced so:
    get_inner_vars(device)[V_pss_var] = K_ω * (ω0 - ω0) + K_p * (ω0 * τe - P0)

end
=#
