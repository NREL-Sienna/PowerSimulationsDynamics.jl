"""
Model of 2-state synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{BaseMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                  A <: AVR,
                                                                                  TG <: TurbineGov,
                                                                                  P <: PSS}



    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]

    #Get parameters
    R = device.machine.R
    Xd_p = device.machine.Xd_p
    eq_p = device.machine.eq_p
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0/(R^2 + Xd_p^2)) * ( Xd_p*(eq_p - V_dq[2]) - R*V_dq[1])  #15.36
    i_q = (1.0/(R^2 + Xd_p^2)) * ( Xd_p*V_dq[1] + R*(eq_p - V_dq[2]) ) #15.36
    Pe = (V_dq[1]+ R*i_d)*i_d + (V_dq[2] + R*i_q)*i_q      #15.35

    #Update inner_vars
    device.inner_vars[τe_var] = Pe #Model assume ω approx 1.0

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end



"""
Model of 4-state (One d- and One q-) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{OneDOneQMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                      A <: AVR,
                                                                                      TG <: TurbineGov,
                                                                                      P <: PSS}




    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]


    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]

    #Get parameters
    R = device.machine.R
    Xd = device.machine.Xd
    Xq = device.machine.Xq
    Xd_p = device.machine.Xd_p
    Xq_p = device.machine.Xq_p
    Td0_p = device.machine.Td0_p
    Tq0_p = device.machine.Tq0_p
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0/(R^2 + Xd_p*Xq_p)) * ( Xq_p*(eq_p - V_dq[2]) + R*(ed_p - V_dq[1]) )  #15.32
    i_q = (1.0/(R^2 + Xd_p*Xq_p)) * ( -Xd_p*(ed_p - V_dq[1]) + R*(eq_p - V_dq[2]) )  #15.32
    Pe = (V_dq[1]+ R*i_d)*i_d + (V_dq[2] + R*i_q)*i_q      #15.35

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0/Td0_p)*(-eq_p - (Xd - Xd_p)*i_d + Vf)     #15.29 eq_p
    output_ode[local_ix[2]] = (1.0/Tq0_p)*(-ed_p + (Xq - Xq_p)*i_q)          #15.30 ed_p

    #Update inner_vars
    device.inner_vars[τe_var] = Pe #We should divide by ω, so we should receive ω

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end


"""
Model of 6-state (MarconatoMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{MarconatoMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                       A <: AVR,
                                                                                       TG <: TurbineGov,
                                                                                       P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]

    #Get parameters
    R = device.machine.R
    Xd = device.machine.Xd
    Xq = device.machine.Xq
    Xd_p = device.machine.Xd_p
    Xq_p = device.machine.Xq_p
    Xd_pp = device.machine.Xd_pp
    Xq_pp = device.machine.Xq_pp
    Td0_p = device.machine.Td0_p
    Tq0_p = device.machine.Tq0_p
    Td0_pp = device.machine.Td0_pp
    Tq0_pp = device.machine.Tq0_pp
    T_AA = device.machine.T_AA
    γd = device.machine.γd
    γq = device.machine.γq
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0/Xd_pp)*(eq_pp - ψd)      #15.18
    i_q = (1.0/Xq_pp)*(-ed_pp - ψq)     #15.18
    τ_e = ψd*i_q - ψq*i_d               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = 2*π*f0*(R*i_q - ω*ψd + V_dq[2])                        #15.9 ψq
    output_ode[local_ix[2]] = 2*π*f0*(R*i_d + ω*ψq + V_dq[1])                        #15.9 ψd
    output_ode[local_ix[3]] = ( (1.0/Td0_p)*(-eq_p - (Xd - Xd_p - γd)*i_d
                                + (1-(T_AA/Td0_p))*Vf) )                             #15.16 eq_p
    output_ode[local_ix[4]] = (1.0/Tq0_p)*(-ed_p + (Xq - Xq_p - γq)*i_q)             #15.16 ed_p
    output_ode[local_ix[5]] = ( (1.0/Td0_pp)*(-eq_pp + eq_p
                                - (Xd_p - Xd_pp + γd)*i_d + (T_AA/Td0_p)*Vf) )       #15.16 eq_pp
    output_ode[local_ix[6]] = (1.0/Tq0_pp)*(-ed_pp + ed_p + (Xq_p - Xq_pp + γq)*i_q) #15.16 ed_pp

    #Update inner_vars
    device.inner_vars[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end



"""
Model of 4-state (SimpleMarconatoMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{SimpleMarconatoMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                             A <: AVR,
                                                                                             TG <: TurbineGov,
                                                                                             P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]

    #Get parameters
    R = device.machine.R
    Xd = device.machine.Xd
    Xq = device.machine.Xq
    Xd_p = device.machine.Xd_p
    Xq_p = device.machine.Xq_p
    Xd_pp = device.machine.Xd_pp
    Xq_pp = device.machine.Xq_pp
    Td0_p = device.machine.Td0_p
    Tq0_p = device.machine.Tq0_p
    Td0_pp = device.machine.Td0_pp
    Tq0_pp = device.machine.Tq0_pp
    T_AA = device.machine.T_AA
    γd = device.machine.γd
    γq = device.machine.γq
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    #i_dq = inv([[R -Xq_pp]; [Xq_pp R]])*[ed_pp - V_dq[1]; eq_pp - V_dq[2]]
    i_d = (1.0/(R^2 + Xd_pp*Xq_pp)) * (Xq_pp*(eq_pp - V_dq[2])  +  R*(ed_pp - V_dq[1]))      #15.25
    i_q = (1.0/(R^2 + Xd_pp*Xq_pp)) * (-Xd_pp*(ed_pp - V_dq[1]) +  R*(eq_pp - V_dq[2]))      #15.25
    τ_e = (V_dq[1]+ R*i_d)*i_d + (V_dq[2] + R*i_q)*i_q               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = ( (1.0/Td0_p)*(-eq_p - (Xd - Xd_p - γd)*i_d
                                + (1-(T_AA/Td0_p))*Vf) )                             #15.16 eq_p
    output_ode[local_ix[2]] = (1.0/Tq0_p)*(-ed_p + (Xq - Xq_p - γq)*i_q)             #15.16 ed_p
    output_ode[local_ix[3]] = ( (1.0/Td0_pp)*(-eq_pp + eq_p
                                - (Xd_p - Xd_pp + γd)*i_d + (T_AA/Td0_p)*Vf) )       #15.16 eq_pp
    output_ode[local_ix[4]] = (1.0/Tq0_pp)*(-ed_pp + ed_p + (Xq_p - Xq_pp + γq)*i_q) #15.16 ed_pp


    #Update inner_vars
    device.inner_vars[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end



"""
Model of 6-state (AndersonFouadMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{AndersonFouadMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                           A <: AVR,
                                                                                           TG <: TurbineGov,
                                                                                           P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]

    #Get parameters
    R = device.machine.R
    Xd = device.machine.Xd
    Xq = device.machine.Xq
    Xd_p = device.machine.Xd_p
    Xq_p = device.machine.Xq_p
    Xd_pp = device.machine.Xd_pp
    Xq_pp = device.machine.Xq_pp
    Td0_p = device.machine.Td0_p
    Tq0_p = device.machine.Tq0_p
    Td0_pp = device.machine.Td0_pp
    Tq0_pp = device.machine.Tq0_pp
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0/Xd_pp)*(eq_pp - ψd)      #15.18
    i_q = (1.0/Xq_pp)*(-ed_pp - ψq)     #15.18
    τ_e = ψd*i_q - ψq*i_d               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = 2*π*f0*(R*i_q - ω*ψd + V_dq[2])                        #15.9 ψq
    output_ode[local_ix[2]] = 2*π*f0*(R*i_d + ω*ψq + V_dq[1])                        #15.9 ψd
    output_ode[local_ix[3]] = (1.0/Td0_p)*(-eq_p + (Xd - Xd_p)*i_d + Vf)             #15.19 eq_p
    output_ode[local_ix[4]] = (1.0/Tq0_p)*(-ed_p + (Xq - Xq_p)*i_q)                  #15.19 ed_p
    output_ode[local_ix[5]] = (1.0/Td0_pp)*(-eq_pp + eq_p - (Xd_p - Xd_pp)*i_d)      #15.19 eq_pp
    output_ode[local_ix[6]] = (1.0/Tq0_pp)*(-ed_pp + ed_p + (Xq_p - Xq_pp)*i_q)      #15.19 ed_pp

    #Update inner_vars
    device.inner_vars[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end




"""
Model of 4-state (SimpleAFMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{SimpleAFMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                           A <: AVR,
                                                                                           TG <: TurbineGov,
                                                                                           P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]

    #Get parameters
    R = device.machine.R
    Xd = device.machine.Xd
    Xq = device.machine.Xq
    Xd_p = device.machine.Xd_p
    Xq_p = device.machine.Xq_p
    Xd_pp = device.machine.Xd_pp
    Xq_pp = device.machine.Xq_pp
    Td0_p = device.machine.Td0_p
    Tq0_p = device.machine.Tq0_p
    Td0_pp = device.machine.Td0_pp
    Tq0_pp = device.machine.Tq0_pp
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0/(R^2 + Xd_pp*Xq_pp)) * (Xq_pp*(eq_pp - V_dq[2])  +  R*(ed_pp - V_dq[1]))      #15.25
    i_q = (1.0/(R^2 + Xd_pp*Xq_pp)) * (-Xd_pp*(ed_pp - V_dq[1]) +  R*(eq_pp - V_dq[2]))      #15.25
    τ_e = (V_dq[1]+ R*i_d)*i_d + (V_dq[2] + R*i_q)*i_q               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0/Td0_p)*(-eq_p + (Xd - Xd_p)*i_d + Vf)            #15.19 eq_p
    output_ode[local_ix[2]] = (1.0/Tq0_p)*(-ed_p + (Xq - Xq_p)*i_q)                 #15.19 ed_p
    output_ode[local_ix[3]] = (1.0/Td0_pp)*(-eq_pp + eq_p - (Xd_p - Xd_pp)*i_d)     #15.19 eq_pp
    output_ode[local_ix[4]] = (1.0/Tq0_pp)*(-ed_pp + ed_p + (Xq_p - Xq_pp)*i_q)     #15.19 ed_pp

    #Update inner_vars
    device.inner_vars[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end



"""
Model of 5-state (FullOrderMachine) synchronous machine in Julia.
Refer to Power System Dynamics: Stability and Control, by J. Machowski, J. Bialek and J. Bumby,
or Power System Stability and Control by P. Kundur, for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{FullMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                           A <: AVR,
                                                                                           TG <: TurbineGov,
                                                                                           P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψd = internal_states[1]
    ψq = internal_states[2]
    ψf = internal_states[1]
    ψ1d = internal_states[2]
    ψ1q = internal_states[3]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]
    #ψd = device.inner_vars[ψd_var]
    #ψq = device.inner_vars[ψq_var]

    #Get parameters
    R = device.machine.R
    R_f = device.machine.R_f
    R_1d= device.machine.R_1d
    R_1q = device.machine.R_1q
    inv_d_fluxlink = device.machine.inv_d_fluxlink
    inv_q_fluxlink = device.machine.inv_q_fluxlink
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric currents] using Flux Linkage equations
    i_df1d = inv_d_fluxlink*[ψd; ψf; ψ1d]  #11.18 in Machowski (3.127, 3.130 and 3.131 in Kundur)
    i_q1q = inv_q_fluxlink*[ψq; ψ1q]        #11.19 in Machowski (3.128 and 3.132 in Kundur)


    #Compute ODEs
    output_ode[local_ix[1]] = 2*pi*f0*(ψq*ω + R*i_df1d[1] + V_dq[1])  #11.32 in Machowski or 3.120 in Kundur
    output_ode[local_ix[2]] = 2*pi*f0*(-ψd*ω + R*i_q1q[1] + V_dq[2])
    output_ode[local_ix[3]] = -R_f*i_df1d[2] + Vf      #11.33 in Machowski or 3.123 in Kundur
    output_ode[local_ix[4]] = -R_1d*i_df1d[3]           #11.33 in Machowski or 3.124 in Kundur
    output_ode[local_ix[5]] = -R_1q*i_q1q[2]            #11.33 in Machowski or 3.125 in Kundur

    #Update stator fluxes and torque
    #ψd = (1.0/ω)*(R*i_q1q[1] + V_dq[2])         #11.32 in Machowski or 3.121 in Kundur
    #ψq = (1.0/ω)*(-R*i_df1d[1] - V_dq[1])       #11.32 in Machowski or 3.120 in Kundur
    τ_e = ψd*i_q1q[1] - ψq*i_df1d[1]            #15.6 in Milano or 3.117 in Kundur

    #Update inner_vars
    device.inner_vars[τe_var] = τ_e
    device.inner_vars[ψd_var] = ψd
    device.inner_vars[ψq_var] = ψq

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_df1d[1]; i_q1q[1]]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end



"""
Model of 3-state (FullOrderMachine) synchronous machine in Julia.
Refer to Power System Dynamics: Stability and Control, by J. Machowski, J. Bialek and J. Bumby,
or Power System Stability and Control by P. Kundur, for the equations
"""

function mdl_machine_ode!(device_states,
                          output_ode,
                          current_r,
                          current_i,
                          Sbase::Float64,
                          f0::Float64,
                          device::DynGenerator{SimpleFullMachine, S, A, TG, P})  where {S <: Shaft,
                                                                                           A <: AVR,
                                                                                           TG <: TurbineGov,
                                                                                           P <: PSS}

    #Obtain indices for component w/r to device
    local_ix = device.local_state_ix[device.machine]

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψf = internal_states[1]
    ψ1d = internal_states[2]
    ψ1q = internal_states[3]

    #Obtain external states inputs for component
    external_ix = device.input_port_mapping[device.machine]
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = device.inner_vars[VR_gen_var]
    V_tI = device.inner_vars[VI_gen_var]
    Vf = device.inner_vars[Vf_var]
    ψd = device.inner_vars[ψd_var]
    ψq = device.inner_vars[ψq_var]

    #Get parameters
    R = device.machine.R
    R_f = device.machine.R_f
    R_1d= device.machine.R_1d
    R_1q = device.machine.R_1q
    inv_d_fluxlink = device.machine.inv_d_fluxlink
    inv_q_fluxlink = device.machine.inv_q_fluxlink
    BaseMVA = device.machine.MVABase

    #RI to dq transformation
    V_dq = ri_dq(δ)*[V_tR; V_tI]

    #Obtain electric currents] using Flux Linkage equations
    i_df1d = inv_d_fluxlink*[ψd; ψf; ψ1d]  #11.18 in Machowski (3.127, 3.130 and 3.131 in Kundur)
    i_q1q = inv_q_fluxlink*[ψq; ψ1q]        #11.19 in Machowski (3.128 and 3.132 in Kundur)


    #Compute ODEs
    output_ode[local_ix[1]] = -R_f*i_df1d[2] + Vf      #11.33 in Machowski or 3.123 in Kundur
    output_ode[local_ix[2]] = -R_1d*i_df1d[3]           #11.33 in Machowski or 3.124 in Kundur
    output_ode[local_ix[3]] = -R_1q*i_q1q[2]            #11.33 in Machowski or 3.125 in Kundur

    #Update stator fluxes and torque
    ψd = (1.0/ω)*(R*i_q1q[1] + V_dq[2])         #11.32 in Machowski or 3.121 in Kundur
    ψq = (1.0/ω)*(-R*i_df1d[1] - V_dq[1])       #11.32 in Machowski or 3.120 in Kundur
    τ_e = ψd*i_q1q[1] - ψq*i_df1d[1]            #15.6 in Milano or 3.117 in Kundur

    #Update inner_vars
    device.inner_vars[τe_var] = τ_e
    device.inner_vars[ψd_var] = ψd
    device.inner_vars[ψq_var] = ψq

    #Compute current from the generator to the grid
    I_RI = (BaseMVA/Sbase)*dq_ri(δ)*[i_df1d[1]; i_q1q[1]]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end
