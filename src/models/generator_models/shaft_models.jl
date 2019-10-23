function mdl_shaft_ode!(device_states,
        output_ode,
        f0::Float64,
        device::DynGenerator{M, SingleMass, A, TG, P})  where {M <: Machine,
                                                               A <: AVR,
                                                               TG <: TurbineGov,
                                                               P <: PSS}



    #Obtain references
    ω_ref = get_ω_ref(device)

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.shaft]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    δ = internal_states[1]
    ω = internal_states[2]

    #Obtain inner variables for component
    τe = device.inner_vars[τe_var]
    τm = device.inner_vars[τm_var]

    #Get parameters
    H = device.shaft.H
    D = device.shaft.D

    #Compute 2 states ODEs
    output_ode[local_ix[1]] = 2*π*f0*(ω-ω_ref)                    #15.5
    output_ode[local_ix[2]] = (1/(2*H))*(τm - τe - D*(ω-ω_ref))   #15.5

    return
end


function mdl_shaft_ode!(device_states,
                        output_ode,
                        f0::Float64,
                        device::DynGenerator{M, FiveMassShaft, A, TG, P})  where {M <: Machine,
                                                                                  A <: AVR,
                                                                                  TG <: TurbineGov,
                                                                                  P <: PSS}

    #Obtain references
    ω_ref = get_ω_ref(device)

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.shaft]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    δ = internal_states[1]
    ω = internal_states[2]
    δ_hp = internal_states[3]
    ω_hp = internal_states[4]
    δ_ip = internal_states[5]
    ω_ip = internal_states[6]
    δ_lp = internal_states[7]
    ω_lp = internal_states[8]
    δ_ex = internal_states[9]
    ω_ex = internal_states[10]

    #Obtain inner variables for component
    τe = device.inner_vars[τe_var]
    τm = device.inner_vars[τm_var]

    #Get parameters
    H = device.shaft.H
    H_hp = device.shaft.H_hp
    H_ip = device.shaft.H_ip
    H_lp = device.shaft.H_lp
    H_ex = device.shaft.H_ex
    D = device.shaft.D
    D_hp = device.shaft.D_hp
    D_ip = device.shaft.D_ip
    D_lp = device.shaft.D_lp
    D_ex = device.shaft.D_ex
    D_12 = device.shaft.D_12
    D_23 = device.shaft.D_23
    D_34 = device.shaft.D_34
    D_45 = device.shaft.D_45
    K_hp = device.shaft.K_hp
    K_ip = device.shaft.K_ip
    K_lp = device.shaft.K_lp
    K_ex = device.shaft.K_ex


    #Compute 10 states ODEs #15.51
    output_ode[local_ix[1]] = 2.0*π*f0*(ω-ω_ref)

    output_ode[local_ix[2]] = (1.0/(2.0*H))*(-τe - D*(ω-ω_ref) - D_34*(ω-ω_lp)
                              -D_45*(ω - ω_ex) + K_lp*(δ_lp - δ)
                              +K_ex*(δ_ex - δ))

    output_ode[local_ix[3]] = 2.0*π*f0*(ω_hp - ω_ref)

    output_ode[local_ix[4]] =  (1.0/(2.0*H_hp))*(τm - D_hp*(ω_hp - ω_ref)
                                - D_12*(ω_hp - ω_ip) +K_hp*(δ_ip - δ_hp))

    output_ode[local_ix[5]] = 2*π*f0*(ω_ip - ω_ref)

    output_ode[local_ix[6]] =  (1.0/(2*H_ip))*(-D_ip*(ω_ip - ω_ref)
                                - D_12*(ω_ip - ω_hp) -D_23*(ω_ip - ω_lp)
                                + K_hp*(δ_hp - δ_ip) +K_ip*(δ_lp - δ_ip))

    output_ode[local_ix[7]] = 2.0*π*f0*(ω_lp - ω_ref)

    output_ode[local_ix[8]] =  (1.0/(2.0*H_lp))*(-D_lp*(ω_lp - ω_ref)
                                - D_23*(ω_lp - ω_ip) -D_34*(ω_lp - ω)
                                + K_ip*(δ_ip - δ_lp) +K_lp*(δ-δ_lp))

    output_ode[local_ix[9]] = 2*π*f0*(ω_ex - ω_ref)

    output_ode[local_ix[10]] = (1.0/(2.0*H_ex))*(-D_ex*(ω_ex - ω_ref)
                                - D_45*(ω_ex - ω) +K_ex*(δ-δ_ex))

    return
end
