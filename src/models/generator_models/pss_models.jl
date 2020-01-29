function mdl_pss_ode!(device_states,
                     output_ode,
                     device::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSFixed})  where {M <: PSY.Machine,
                                                                        S <: PSY.Shaft,
                                                                        A <: PSY.AVR,
                                                                        TG <: PSY.TurbineGov}

    #Update V_pss on inner vars
    get_inner_vars(device)[V_pss_var] = PSY.get_V_pss(PSY.get_pss(device))

    return
end



function mdl_pss_ode!(device_states,
                     output_ode,
                     device::PSY.DynamicGenerator{M, S, A, TG, PSY.PSSSimple})  where {M <: PSY.Machine,
                                                                        S <: PSY.Shaft,
                                                                        A <: PSY.AVR,
                                                                        TG <: PSY.TurbineGov}

    #Get references
    ω_ref = PSY.get_ω_ref(device)
    P_ref = PSY.get_P_ref(device)

    #Obtain external states for device
    external_ix = get_input_port_ix(device, PSY.PSSSimple)
    ω = @view device_states[external_ix]

    #Define external inner vars for component
    V_th = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)
    τe = get_inner_vars(device)[τe_var]

    #Get parameters
    pss = PSY.get_pss(device)
    K_ω = PSY.get_K_ω(pss)
    K_p = PSY.get_K_p(pss)

    #Update V_pss on inner vars
    get_inner_vars(device)[V_pss_var] = K_ω*(ω-ω_ref) + K_p*(ω*τe - P_ref)

    return
end
