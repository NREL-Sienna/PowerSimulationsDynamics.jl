function mdl_pss_ode!(
    device_states,
    output_ode,
    ω_sys,
    dyn_data::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Update V_pss on inner vars
    get_inner_vars(dyn_data)[V_pss_var] = PSY.get_V_pss(PSY.get_pss(dyn_data))

    return
end

#Currently not working properly.
#=
function mdl_pss_ode!(
    device_states,
    output_ode,
    ω_sys,
    dyn_data::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple},
) where {M <: PSY.Machine, S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov}

    #Get references
    P0 = PSY.get_active_power(static)

    #Obtain external states for device
    external_ix = get_input_port_ix(dyn_data, PSY.PSSSimple)
    ω = @view device_states[external_ix]

    #Define external inner vars for component
    V_th = sqrt(get_inner_vars(dyn_data)[VR_gen_var]^2 + get_inner_vars(dyn_data)[VI_gen_var]^2)
    τe = get_inner_vars(dyn_data)[τe_var]

    #Get parameters
    pss = PSY.get_pss(dyn_data)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #Update V_pss on inner vars
    get_inner_vars(dyn_data)[V_pss_var] = K_ω * (ω[1] - ω_sys) + K_p * (ω[1] * τe - P0)

    return
end
=#
