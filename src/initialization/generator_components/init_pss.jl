function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dyn_data::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov} end

#Currently not working properly.
#=
function initialize_pss!(
    device_states,
    static::PSY.StaticInjection,
    dyn_data::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    ω0 = PSY.get_ext(dyn_data)[CONTROL_REFS][ω_ref_index]
    P0 = PSY.get_active_power(static)
    τe = get_inner_vars(dyn_data)[τe_var]

    #Get parameters
    pss = PSY.get_pss(dyn_data)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #It assumes that the system is balanced so:
    get_inner_vars(dyn_data)[V_pss_var] = K_ω * (ω0 - ω0) + K_p * (ω0 * τe - P0)

end
=#
