"""
Model of 0-state synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.BaseMachine)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    eq_p = PSY.get_eq_p(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (eq_p - V_dq[2]) - R * V_dq[1])  #15.36
    i_q = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_dq[1] + R * (eq_p - V_dq[2])) #15.36
    Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q      #15.35

    #Update inner_vars
    get_inner_vars(device)[τe_var] = Pe #Model assume ω approx 1.0

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

"""
Model of 2-state (One d- and One q-) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.OneDOneQMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.OneDOneQMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.OneDOneQMachine)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0 / (R^2 + Xd_p * Xq_p)) * (Xq_p * (eq_p - V_dq[2]) + R * (ed_p - V_dq[1]))  #15.32
    i_q = (1.0 / (R^2 + Xd_p * Xq_p)) * (-Xd_p * (ed_p - V_dq[1]) + R * (eq_p - V_dq[2]))  #15.32
    Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q      #15.35

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (-eq_p - (Xd - Xd_p) * i_d + Vf)     #15.29 eq_p
    output_ode[local_ix[2]] = (1.0 / Tq0_p) * (-ed_p + (Xq - Xq_p) * i_q)          #15.30 ed_p

    #Update inner_vars
    get_inner_vars(device)[τe_var] = Pe #Model assume ω approx 1.0

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

"""
Model of 6-state (MarconatoMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.MarconatoMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.MarconatoMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.MarconatoMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    T_AA = PSY.get_T_AA(machine)
    γd = PSY.get_γd(machine)
    γq = PSY.get_γq(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0 / Xd_pp) * (eq_pp - ψd)      #15.18
    i_q = (1.0 / Xq_pp) * (-ed_pp - ψq)     #15.18
    τ_e = ψd * i_q - ψq * i_d               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = 2 * π * f0 * (R * i_q - ω * ψd + V_dq[2])                        #15.9 ψq
    output_ode[local_ix[2]] = 2 * π * f0 * (R * i_d + ω * ψq + V_dq[1])                        #15.9 ψd
    output_ode[local_ix[3]] =
        ((1.0 / Td0_p) * (-eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf))                             #15.16 eq_p
    output_ode[local_ix[4]] = (1.0 / Tq0_p) * (-ed_p + (Xq - Xq_p - γq) * i_q)             #15.16 ed_p
    output_ode[local_ix[5]] =
        ((1.0 / Td0_pp) * (-eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf))       #15.16 eq_pp
    output_ode[local_ix[6]] = (1.0 / Tq0_pp) * (-ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q) #15.16 ed_pp

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

"""
Model of 4-state (SimpleMarconatoMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.SimpleMarconatoMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.SimpleMarconatoMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.SimpleMarconatoMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    T_AA = PSY.get_T_AA(machine)
    γd = PSY.get_γd(machine)
    γq = PSY.get_γq(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    #i_dq = inv([[R -Xq_pp]; [Xq_pp R]])*[ed_pp - V_dq[1]; eq_pp - V_dq[2]]
    i_d =
        (1.0 / (R^2 + Xd_pp * Xq_pp)) * (Xq_pp * (eq_pp - V_dq[2]) + R * (ed_pp - V_dq[1]))      #15.25
    i_q =
        (1.0 / (R^2 + Xd_pp * Xq_pp)) * (-Xd_pp * (ed_pp - V_dq[1]) + R * (eq_pp - V_dq[2]))      #15.25
    τ_e = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] =
        ((1.0 / Td0_p) * (-eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf))                             #15.16 eq_p
    output_ode[local_ix[2]] = (1.0 / Tq0_p) * (-ed_p + (Xq - Xq_p - γq) * i_q)             #15.16 ed_p
    output_ode[local_ix[3]] =
        ((1.0 / Td0_pp) * (-eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf))       #15.16 eq_pp
    output_ode[local_ix[4]] = (1.0 / Tq0_pp) * (-ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q) #15.16 ed_pp

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

"""
Model of 6-state (AndersonFouadMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.AndersonFouadMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.AndersonFouadMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.AndersonFouadMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0 / Xd_pp) * (eq_pp - ψd)      #15.18
    i_q = (1.0 / Xq_pp) * (-ed_pp - ψq)     #15.18
    τ_e = ψd * i_q - ψq * i_d               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = 2 * π * f0 * (R * i_q - ω * ψd + V_dq[2])                        #15.9 ψq
    output_ode[local_ix[2]] = 2 * π * f0 * (R * i_d + ω * ψq + V_dq[1])                        #15.9 ψd
    output_ode[local_ix[3]] = (1.0 / Td0_p) * (-eq_p + (Xd - Xd_p) * i_d + Vf)             #15.19 eq_p
    output_ode[local_ix[4]] = (1.0 / Tq0_p) * (-ed_p + (Xq - Xq_p) * i_q)                  #15.19 ed_p
    output_ode[local_ix[5]] = (1.0 / Td0_pp) * (-eq_pp + eq_p - (Xd_p - Xd_pp) * i_d)      #15.19 eq_pp
    output_ode[local_ix[6]] = (1.0 / Tq0_pp) * (-ed_pp + ed_p + (Xq_p - Xq_pp) * i_q)      #15.19 ed_pp

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

"""
Model of 4-state (SimpleAFMachine) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.SimpleAFMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.SimpleAFMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.SimpleAFMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d =
        (1.0 / (R^2 + Xd_pp * Xq_pp)) * (Xq_pp * (eq_pp - V_dq[2]) + R * (ed_pp - V_dq[1]))      #15.25
    i_q =
        (1.0 / (R^2 + Xd_pp * Xq_pp)) * (-Xd_pp * (ed_pp - V_dq[1]) + R * (eq_pp - V_dq[2]))      #15.25
    τ_e = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q               #15.6

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (-eq_p + (Xd - Xd_p) * i_d + Vf)            #15.19 eq_p
    output_ode[local_ix[2]] = (1.0 / Tq0_p) * (-ed_p + (Xq - Xq_p) * i_q)                 #15.19 ed_p
    output_ode[local_ix[3]] = (1.0 / Td0_pp) * (-eq_pp + eq_p - (Xd_p - Xd_pp) * i_d)     #15.19 eq_pp
    output_ode[local_ix[4]] = (1.0 / Tq0_pp) * (-ed_pp + ed_p + (Xq_p - Xq_pp) * i_q)     #15.19 ed_pp

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end


"""
Model of 4-state (RoundRotorQuadratic - GENROU) synchronous machine in Julia.
Refer to SynchGen and Excitation Models by Paszek et al. for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.RoundRotorQuadratic, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.RoundRotorQuadratic)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    ψ_kd = internal_states[3]
    ψ_kq = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.RoundRotorQuadratic)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var] #E_fd: Field voltage

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    Asat, Bsat = PSY.get_saturation_coeffs(machine)
    basepower = PSY.get_base_power(device)

    #Additional Parameters
    γ_d1 = (Xd_pp - Xl) / (Xd_p - Xl)
    γ_q1 = (Xq_pp - Xl) / (Xq_p - Xl)
    γ_d2 = (Xd_p - Xd_pp) / (Xd_p - Xl)^2
    γ_q2 = (Xq_p - Xq_pp) / (Xq_p - Xl)^2

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]
    
    #Obtain electric variables
    i_d = (1.0 / (R^2 + Xd_pp * Xq_pp) ) * ( R * (ed_p * γ_q1 - ψ_kq * (1 - γ_q1) - V_dq[1]) - Xq_pp * (V_dq[2] - eq_p * γ_d1 - ψ_kd * (1 - γ_d1)) ) #2.21a
    i_q = (1.0 / (R^2 + Xd_pp * Xq_pp) ) * ( - R * (V_dq[2] - eq_p * γ_d1 - ψ_kd * (1 - γ_d1)) -  Xd_pp * (ed_p * γ_q1 - V_dq[1] - ψ_kq * ( 1 - γ_d1)) ) #2.20a
    τ_e = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q         #2.36c + 2.34a + 2.35a
    
    #Get auxiliary variables
    #dtψ_kd = (1.0 / Td0_pp) * (eq_p - ψ_kd - i_d*(Xd_p - Xl)) #2.20c dψ_kd/dt
    #dtψ_kq = -(1.0 / Tq0_pp) * (-ed_p + ψ_kq + i_q*(Xq_p - Xl))
    #ΔId = Se_ψ * ψd_pp / ψ_pp
    #ΔIq = Se_ψ * (ψq_pp / ψ_pp) * ((Xq - Xl)/(Xd - Xl))
    #Xad_Ifd = eq_p + (Xd - Xd_p)*(i_d + Td0_pp * dtψ_kd * (Xd_p - Xd_pp)/(Xd_p - Xl)^2) + ΔId
    #Xaq_Ikq = ed_p + (Xq - Xq_p')*(i_q + Tq0_pp * dtψ_kq * (Xq_p - Xq_pp)/(Xq_p - Xl)^2) + ΔIq

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - eq_p - (Xd - Xd_p) * (i_d - γ_d2 * ψ_kd - (1 - γ_d1)*i_d + γ_d2 * eq_p) )                        #15.13 eq_p
    output_ode[local_ix[2]] = (1.0 / Td0_p) * (- ed_p + (Xq - Xq_p) * (i_q - γ_q2 * ψ_kq - (1 - γ_q1)*i_q - γ_d2 * ed_p ) )                    #15.13 ed_p
    output_ode[local_ix[3]] = (1.0 / Td0_pp) * (eq_p - ψ_kd - (Xd_p - Xl) * i_d)                    #15.13 ψ_kd
    output_ode[local_ix[4]] = (1.0 / Tq0_pp) * (-ed_p - ψ_kq - (Xq_p - Xl) * i_q)                 #15.19 ψ_kq

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end


#=
"""
Model of 4-state (RoundRotorQuadratic - GENROU) synchronous machine in Julia.
Refer to SynchGen and Excitation Models by Paszek et al. for the equations
"""
function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.RoundRotorQuadratic, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.RoundRotorQuadratic)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    ψ_kd = internal_states[3]
    ψ_kq = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.RoundRotorQuadratic)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var] #E_fd: Field voltage

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    Asat, Bsat = PSY.get_saturation_coeffs(machine)
    basepower = PSY.get_base_power(device)

    #Additional Parameters
    γ_d1 = (Xd_pp - Xl) / (Xd_p - Xl)
    γ_q1 = (Xq_pp - Xl) / (Xq_p - Xl)
    γ_d2 = (Xd_p - Xd_pp) / (Xd_p - Xl)^2
    γ_q2 = (Xq_p - Xq_pp) / (Xq_p - Xl)^2

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]
    
    #Compute additional terms
    ψd_pp = eq_p * (Xd_pp - Xl)/(Xd_p - Xl) + ψ_kd * (Xd_p - Xd_pp)/(Xd_p - Xl) #2.20d
    ψq_pp = ed_p * (Xq_pp - Xl)/(Xq_p - Xl) + ψ_kq * (Xq_p - Xq_pp)/(Xq_p - Xl) #2.21d
    ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)
    Se_ψ = Bsat * (ψ_pp - Asat)^2

    #Obtain electric variables
    i_d = (1.0 / (R^2 + Xd_pp * Xq_pp) ) * ( Xq_pp * (ψd_pp - V_dq[2]) - R * (V_dq[1] + ψq_pp)) #2.21a
    i_q = (1.0 / (R^2 + Xd_pp * Xq_pp) ) * ( Xd_pp * (ψq_pp + V_dq[1]) - R * (V_dq[2] - ψd_pp)) #2.20a
    τ_e = ψd_pp * i_q - ψq_pp * i_d         #2.36c + 2.34a + 2.35a
    
    #Get auxiliary variables
    dtψ_kd = (1.0 / Td0_pp) * (eq_p - ψ_kd - i_d*(Xd_p - Xl)) #2.20c dψ_kd/dt
    dtψ_kq = -(1.0 / Tq0_pp) * (-ed_p + ψ_kq + i_q*(Xq_p - Xl))
    ΔId = Se_ψ * ψd_pp / ψ_pp
    ΔIq = Se_ψ * (ψq_pp / ψ_pp) * ((Xq - Xl)/(Xd - Xl))
    Xad_Ifd = eq_p + (Xd - Xd_p)*(i_d + Td0_pp * dtψ_kd * (Xd_p - Xd_pp)/(Xd_p - Xl)^2) + ΔId
    Xaq_Ikq = ed_p + (Xq - Xq_p')*(i_q + Tq0_pp * dtψ_kq * (Xq_p - Xq_pp)/(Xq_p - Xl)^2) + ΔIq

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - Xad_Ifd)                        #2.20 eq_p
    output_ode[local_ix[2]] = (1.0 / Td0_p) * (- Xaq_Ikq)                    #15.9 ed_p
    output_ode[local_ix[3]] = dtψ_kd                     #2.20c ψ_kd
    output_ode[local_ix[4]] = dtψ_kq           #15.19 ψ_kq

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end
=#


#Not already implemented:
#=
"""
Model of 5-state (FullOrderMachine) synchronous machine in Julia.
Refer to Power System Dynamics: Stability and Control, by J. Machowski, J. Bialek and J. Bumby,
or Power System Stability and Control by P. Kundur, for the equations
"""

function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.FullMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.FullMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψd = internal_states[1]
    ψq = internal_states[2]
    ψf = internal_states[3]
    ψ1d = internal_states[4]
    ψ1q = internal_states[5]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.FullMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]
    #ψd = get_inner_vars(device)[ψd_var]
    #ψq = get_inner_vars(device)[ψq_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    R_f = PSY.get_R_f(machine)
    R_1d = PSY.get_R_1d(machine)
    R_1q = PSY.get_R_1q(machine)
    inv_d_fluxlink = PSY.get_inv_d_fluxlink(machine)
    inv_q_fluxlink = PSY.get_inv_q_fluxlink(machine)
    basepower = PSY.get_base_power(device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric currents] using Flux Linkage equations
    i_df1d = inv_d_fluxlink * [ψd; ψf; ψ1d]  #11.18 in Machowski (3.127, 3.130 and 3.131 in Kundur)
    i_q1q = inv_q_fluxlink * [ψq; ψ1q]        #11.19 in Machowski (3.128 and 3.132 in Kundur)

    #Compute ODEs
    output_ode[local_ix[1]] = 2 * pi * f0 * (ψq * ω + R * i_df1d[1] + V_dq[1])  #11.32 in Machowski or 3.120 in Kundur
    output_ode[local_ix[2]] = 2 * pi * f0 * (-ψd * ω + R * i_q1q[1] + V_dq[2])
    output_ode[local_ix[3]] = -R_f * i_df1d[2] + Vf      #11.33 in Machowski or 3.123 in Kundur
    output_ode[local_ix[4]] = -R_1d * i_df1d[3]           #11.33 in Machowski or 3.124 in Kundur
    output_ode[local_ix[5]] = -R_1q * i_q1q[2]            #11.33 in Machowski or 3.125 in Kundur

    #Update stator fluxes and torque
    #ψd = (1.0/ω)*(R*i_q1q[1] + V_dq[2])         #11.32 in Machowski or 3.121 in Kundur
    #ψq = (1.0/ω)*(-R*i_df1d[1] - V_dq[1])       #11.32 in Machowski or 3.120 in Kundur
    τ_e = ψd * i_q1q[1] - ψq * i_df1d[1]            #15.6 in Milano or 3.117 in Kundur

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e
    get_inner_vars(device)[ψd_var] = ψd
    get_inner_vars(device)[ψq_var] = ψq

    #Compute current from the generator to the grid
    I_RI = basepower * dq_ri(δ) * [i_df1d[1]; i_q1q[1]]

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

function mdl_machine_ode!(
    device_states,
    output_ode,
    current_r,
    current_i,
    Sbase::Float64,
    f0::Float64,
    device::PSY.DynamicGenerator{PSY.SimpleFullMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.SimpleFullMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψf = internal_states[1]
    ψ1d = internal_states[2]
    ψ1q = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(device, PSY.SimpleFullMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = get_inner_vars(device)[VR_gen_var]
    V_tI = get_inner_vars(device)[VI_gen_var]
    Vf = get_inner_vars(device)[Vf_var]
    ψd = get_inner_vars(device)[ψd_var]
    ψq = get_inner_vars(device)[ψq_var]

    #Get parameters
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    R_f = PSY.get_R_f(machine)
    R_1d = PSY.get_R_1d(machine)
    R_1q = PSY.get_R_1q(machine)
    inv_d_fluxlink = PSY.get_inv_d_fluxlink(machine)
    inv_q_fluxlink = PSY.get_inv_q_fluxlink(machine)
    BaseMVA = PSY.get_MVABase(machine)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric currents] using Flux Linkage equations
    i_df1d = inv_d_fluxlink * [ψd; ψf; ψ1d]  #11.18 in Machowski (3.127, 3.130 and 3.131 in Kundur)
    i_q1q = inv_q_fluxlink * [ψq; ψ1q]        #11.19 in Machowski (3.128 and 3.132 in Kundur)

    #Compute ODEs
    output_ode[local_ix[1]] = -R_f * i_df1d[2] + Vf      #11.33 in Machowski or 3.123 in Kundur
    output_ode[local_ix[2]] = -R_1d * i_df1d[3]           #11.33 in Machowski or 3.124 in Kundur
    output_ode[local_ix[3]] = -R_1q * i_q1q[2]            #11.33 in Machowski or 3.125 in Kundur

    #Update stator fluxes and torque
    ψd = (1.0 / ω) * (R * i_q1q[1] + V_dq[2])         #11.32 in Machowski or 3.121 in Kundur
    ψq = (1.0 / ω) * (-R * i_df1d[1] - V_dq[1])       #11.32 in Machowski or 3.120 in Kundur
    τ_e = ψd * i_q1q[1] - ψq * i_df1d[1]            #15.6 in Milano or 3.117 in Kundur

    #Update inner_vars
    get_inner_vars(device)[τe_var] = τ_e
    get_inner_vars(device)[ψd_var] = ψd
    get_inner_vars(device)[ψq_var] = ψq

    #Compute current from the generator to the grid
    I_RI = (BaseMVA / Sbase) * dq_ri(δ) * [i_df1d[1]; i_q1q[1]]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end
=#
