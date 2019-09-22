function mdl_tg_ode!(device_states,
                    output_ode,
                    device::DynGenerator{M, S, A, TGFixed, P})  where {M <: Machine,
                                                                       S <: Shaft,
                                                                       A <: AVR,
                                                                       P <: PSS}

    #Update inner vars
    device.inner_vars[τm_var] = get_P_ref(device) * device.tg.efficiency
    #@show device
    #@show device.inner_vars[τm_var]


    return

end


function mdl_tg_ode!(device_states,
                    output_ode,
                    device::DynGenerator{M, S, A, TGTypeI, P})  where {M <: Machine,
                                                                       S <: Shaft,
                                                                       A <: AVR,
                                                                       P <: PSS}

    #Obtain references
    ω_ref = get_ω_ref(device)
    P_ref = get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.tg]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    x_g1 = internal_states[1]
    x_g2 = internal_states[2]
    x_g3 = internal_states[3]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.tg]
    ω = @view device_states[external_ix]

    #Get Parameters
    R = device.tg.R
    Ts = device.tg.Ts
    Tc = device.tg.Tc
    T3 = device.tg.T3
    T4 = device.tg.T4
    T5 = device.tg.T5
    P_min = device.tg.P_min
    P_max = device.tg.P_max

    #Compute auxiliary parameters
    P_in = P_ref + (1.0/R)*(ω_ref - ω[1])

    #Set anti-windup for P_in. NOT WORKING
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
    device.inner_vars[τm_var] = x_g3 + (T4/T5)*(x_g2 + (T3/Tc)*x_g1)

    return
end




function mdl_tg_ode!(device_states,
                    output_ode,
                    device::DynGenerator{M, S, A, TGTypeII, P})  where {M <: Machine,
                                                                       S <: Shaft,
                                                                       A <: AVR,
                                                                       P <: PSS}

    #Obtain references
    ω_ref = get_ω_ref(device)
    P_ref = get_P_ref(device)

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.tg]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    xg = internal_states[1]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.tg]
    ω = @view device_states[external_ix]

    #Get Parameters
    R = device.tg.R
    T1 = device.tg.T1
    T2 = device.tg.T2
    τ_min = device.tg.τ_min
    τ_max = device.tg.τ_max

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
    device.inner_vars[τm_var] = τ_m

    return
end
