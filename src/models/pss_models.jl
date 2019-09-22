function mdl_pss_ode!(device_states,
                     output_ode,
                     device::DynGenerator{M, S, A, TG, PSSFixed})  where {M <: Machine,
                                                                        S <: Shaft,
                                                                        A <: AVR,
                                                                        TG <: TurbineGov}

    #Compute Vs update
    device.inner_vars[V_pss_var] = device.pss.Vs

    return
end



function mdl_pss_ode!(device_states,
                     output_ode,
                     device::DynGenerator{M, S, A, TG, PSSSimple})  where {M <: Machine,
                                                                        S <: Shaft,
                                                                        A <: AVR,
                                                                        TG <: TurbineGov}

    #Get references
    ω_ref = get_ω_ref(device)
    P_ref = get_P_ref(device)

    #Obtain external states for device
    external_ix = device.input_port_mapping[device.tg]
    ω = @view device_states[external_ix]

    #Define external inner vars for component
    V_th = device.inner_vars[Vh_var]
    τe = device.inner_vars[τe_var]

    #Get parameters
    K_ω = device.pss.K_ω
    K_p = device.pss.K_p

    #Update inner vars
    device.inner_vars[V_pss_var] = K_ω*(ω-ω_ref) + K_p*(ω*τe - P_ref)




    return
end
