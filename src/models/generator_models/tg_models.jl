function mdl_tg_ode!(device_states,
                    output_ode,
                    device::PSY.DynamicGenerator{M, S, A, PSY.TGFixed, P})  where {M <: PSY.Machine,
                                                                       S <: PSY.Shaft,
                                                                       A <: PSY.AVR,
                                                                       P <: PSY.PSS}

    #Update inner vars
    get_inner_vars(device)[τm_var] = PSY.get_P_ref(device) * PSY.get_efficiency(PSY.get_tg(device))

    return
end


function mdl_tg_ode!(device_states,
                    output_ode,
                    device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeI, P})  where {M <: PSY.Machine,
                                                                       S <: PSY.Shaft,
                                                                       A <: PSY.AVR,
                                                                       P <: PSY.PSS}

    #Obtain references
    ω_ref = PSY.get_ω_ref(device)
    P_ref = PSY.get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.TGTypeI)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    x_g3 = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.TGTypeI)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_tg(device)
    R = PSY.get_R(tg)
    Ts = PSY.get_Ts(tg)
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    P_min = PSY.get_P_min(tg)
    P_max = PSY.get_P_max(tg)

    #Compute auxiliary parameters
    P_in = P_ref + (1.0/R)*(ω_ref - ω[1])

    #Set anti-windup for P_in. #TODO in callbacks
    #if P_in > P_max
    #    P_in = P_max
    #elseif P_in < P_min
    #    P_in = P_min
    #end

    #Compute 3 States TG ODE:
    output_ode[local_ix[1]] = (1.0/Ts)*(P_in - x_g1)
    output_ode[local_ix[2]] = (1.0/Tc)*((1.0 - T3/Tc)*x_g1 - x_g2)
    output_ode[local_ix[3]] = (1.0/T5)*((1.0 - T4/T5)*(x_g2 + (T3/Tc)*x_g1) - x_g3)

    #Update mechanical torque
    get_inner_vars(device)[τm_var] = x_g3 + (T4/T5)*(x_g2 + (T3/Tc)*x_g1)

    return
end




function mdl_tg_ode!(device_states,
                    output_ode,
                    device::PSY.DynamicGenerator{M, S, A, PSY.TGTypeII, P})  where {M <: PSY.Machine,
                                                                       S <: PSY.Shaft,
                                                                       A <: PSY.AVR,
                                                                       P <: PSY.PSS}

    #Obtain references
    ω_ref = PSY.get_ω_ref(device)
    P_ref = PSY.get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.TGTypeII)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    xg = internal_states[1]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.TGTypeII)
    ω = @view device_states[external_ix]

    #Get Parameters
    tg = PSY.get_tg(device)
    R = PSY.get_R(tg)
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    τ_min = PSY.get_τ_min(tg)
    τ_max = PSY.get_τ_max(tg)

    #Compute auxiliary parameters
    τ_m = xg + (1.0/R)*(ω_ref - ω[1]) + P_ref/1.0

    #Set anti-windup for τ_m. NOT WORKING
    #if τ_m > τ_max
    #    τ_m = τ_max
    #elseif τ_m < τ_min
    #    τ_m = τ_min
    #end

    #Compute 1 State TG ODE:
    output_ode[local_ix[1]] = (1.0/T2)*( (1.0/R)*(1- T2/T1)*(ω_ref - ω[1]) - xg )

    #Update mechanical torque
    get_inner_vars(device)[τm_var] = τ_m

    return
end
