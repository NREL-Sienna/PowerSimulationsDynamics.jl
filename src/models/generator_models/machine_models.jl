function mass_matrix_machine_entries!(
    mass_matrix,
    machine::M,
    global_index::Base.ImmutableDict{Symbol, Int64},
) where {M <: PSY.Machine}
    @debug "Using default mass matrix entries $M"
end

"""
Model of 0-state synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function mdl_machine_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.BaseMachine)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    eq_p = PSY.get_eq_p(machine)
    basepower = PSY.get_base_power(dynamic_device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Obtain electric variables
    i_d = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (eq_p - V_dq[2]) - R * V_dq[1])  #15.36
    i_q = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_dq[1] + R * (eq_p - V_dq[2])) #15.36
    Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q      #15.35

    #Update inner_vars
    inner_vars[τe_var] = Pe #Model assume ω approx 1.0

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
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.OneDOneQMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.OneDOneQMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.OneDOneQMachine)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Tq0_p = PSY.get_Tq0_p(machine)
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = Pe #Model assume ω approx 1.0

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
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.MarconatoMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.MarconatoMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.MarconatoMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = τ_e

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
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SimpleMarconatoMachine, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)
    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SimpleMarconatoMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.SimpleMarconatoMachine)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = τ_e

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
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.AndersonFouadMachine, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.AndersonFouadMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψq = internal_states[1]
    ψd = internal_states[2]
    eq_p = internal_states[3]
    ed_p = internal_states[4]
    eq_pp = internal_states[5]
    ed_pp = internal_states[6]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.AndersonFouadMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = τ_e

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
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.SimpleAFMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SimpleAFMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    eq_pp = internal_states[3]
    ed_pp = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.SimpleAFMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = τ_e

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [i_d; i_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

function mdl_machine_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
) where {
    M <: Union{PSY.RoundRotorQuadratic, PSY.RoundRotorExponential},
    S <: PSY.Shaft,
    A <: PSY.AVR,
    TG <: PSY.TurbineGov,
    P <: PSY.PSS,
}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, M)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    ψ_kd = internal_states[3]
    ψ_kq = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, M)
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var] #E_fd: Field voltage

    #Get parameters
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
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    γ_q2 = PSY.get_γ_q2(machine)
    γ_qd = PSY.get_γ_qd(machine)
    basepower = PSY.get_base_power(dynamic_device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψq_pp = γ_q1 * ed_p + ψ_kq * (1 - γ_q1)
    ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψ_kd
    ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)
    #Currents
    I_d =
        (1.0 / (R^2 + Xq_pp * Xd_pp)) *
        (-R * (V_dq[1] - ψq_pp) + Xq_pp * (-V_dq[2] + ψd_pp))
    I_q =
        (1.0 / (R^2 + Xq_pp * Xd_pp)) * (Xd_pp * (V_dq[1] - ψq_pp) + R * (-V_dq[2] + ψd_pp))
    Se = saturation_function(machine, ψ_pp)
    Xad_Ifd = eq_p + (Xd - Xd_p) * (γ_d1 * I_d - γ_d2 * ψ_kd + γ_d2 * eq_p) + Se * ψd_pp
    Xaq_I1q =
        ed_p + (Xq - Xq_p) * (γ_q2 * ed_p - γ_q2 * ψ_kq - γ_q1 * I_q) + Se * ψq_pp * γ_qd
    τ_e = I_d * (V_dq[1] + I_d * R) + I_q * (V_dq[2] + I_q * R)

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - Xad_Ifd)                        #2.20b eq_p
    output_ode[local_ix[2]] = (1.0 / Tq0_p) * (-Xaq_I1q)                            #2.21b ed_p
    output_ode[local_ix[3]] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d)   #2.20c ψ_kd
    output_ode[local_ix[4]] = (1.0 / Tq0_pp) * (-ψ_kq + ed_p + (Xq_p - Xl) * I_q)   #2.21c ψ_kq

    #Update inner_vars
    inner_vars[τe_var] = τ_e
    inner_vars[Xad_Ifd_var] = Xad_Ifd

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [I_d; I_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

function mdl_machine_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleQuadratic, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, typeof(machine))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ψ_kd = internal_states[2]
    ψq_pp = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(machine))
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var] #E_fd: Field voltage

    #Get parameters
    R = PSY.get_R(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    basepower = PSY.get_base_power(dynamic_device)

    #RI to dq transformation
    V_d, V_q = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd

    #Currents
    I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d + ψq_pp) + Xd_pp * (ψd_pp - V_q))
    I_q = (1.0 / (R^2 + Xd_pp^2)) * (Xd_pp * (V_d + ψq_pp) + R * (ψd_pp - V_q))
    Se = saturation_function(machine, eq_p)
    Xad_Ifd =
        eq_p + Se * eq_p + (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))
    τ_e = I_d * (V_d + I_d * R) + I_q * (V_q + I_q * R)

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - Xad_Ifd)                        #2.34b eq_p
    output_ode[local_ix[2]] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d)   #2.34c ψ_kd
    output_ode[local_ix[3]] = (1.0 / Tq0_pp) * (-ψq_pp - (Xq - Xq_pp) * I_q)        #2.35b ψq_pp

    #Update inner_vars
    inner_vars[τe_var] = τ_e
    inner_vars[Xad_Ifd_var] = Xad_Ifd

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [I_d; I_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

function mdl_machine_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleExponential, S, A, TG, P},
    },
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, typeof(machine))

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ψ_kd = internal_states[2]
    ψq_pp = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, typeof(machine))
    δ = device_states[external_ix[1]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var] #E_fd: Field voltage

    #Get parameters
    R = PSY.get_R(machine)
    Td0_p = PSY.get_Td0_p(machine)
    Td0_pp = PSY.get_Td0_pp(machine)
    Tq0_pp = PSY.get_Tq0_pp(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    γ_qd = (Xq - Xl) / (Xd - Xl)
    basepower = PSY.get_base_power(dynamic_device)

    #RI to dq transformation
    V_d, V_q = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd
    ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)

    #Currents
    I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d - ψq_pp) + Xq_pp * (-V_q + ψd_pp))
    I_q = (1.0 / (R^2 + Xd_pp^2)) * (Xd_pp * (V_d - ψq_pp) + R * (-V_q + ψd_pp))
    Se = saturation_function(machine, ψ_pp)
    Xad_Ifd =
        eq_p + Se * ψd_pp + (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))
    τ_e = I_d * (V_d + I_d * R) + I_q * (V_q + I_q * R)

    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - Xad_Ifd)                        #2.34b eq_p
    output_ode[local_ix[2]] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d)   #2.34c ψ_kd
    output_ode[local_ix[3]] =
        (1.0 / Tq0_pp) * (-ψq_pp + (Xq - Xq_pp) * I_q - Se * γ_qd * ψq_pp)        #2.35b ψq_pp

    #Update inner_vars
    inner_vars[τe_var] = τ_e
    inner_vars[Xad_Ifd_var] = Xad_Ifd

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [I_d; I_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

function mdl_machine_ode!(
    device_states::AbstractArray{<:ACCEPTED_REAL_TYPES},
    output_ode::AbstractArray{<:ACCEPTED_REAL_TYPES},
    inner_vars::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.GENQEC, S, A, TG, P}},
    ) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    Sbase = get_system_base_power(dynamic_device)
    f0 = get_system_base_frequency(dynamic_device)

    #Obtain indices for component w/r to device
    machine = PSY.get_machine(dynamic_device)
    local_ix = get_local_state_ix(dynamic_device, M)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    eq_p = internal_states[1]
    ed_p = internal_states[2]
    ψd_p = internal_states[3]
    ψq_p = internal_states[4]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, M)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]] 

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var] #E_fd: Field voltage

    #Get parameters
    sat_flag = PSY.get_sat_flag(machine)
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
    Xq_pp = PSY.get_Xq_pp(machine)
    Xl = PSY.get_Xl(machine)
    Kw = PSY.get_Kw(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    γ_q2 = PSY.get_γ_q2(machine)
    basepower = PSY.get_base_power(dynamic_device)

    #RI to dq transformation
    V_dq = ri_dq(δ) * [V_tR; V_tI]

    #Additional Fluxes
    
    ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψd_p 
    ψq_pp = - γ_q1 * ed_p - γ_q2 * (Xq_p - Xl) * ψq_p 
    
    # Get available currents and convert to machine base
    I_R = current_r[1] * (Sbase / basepower)
    I_I = current_i[1] * (Sbase / basepower)
    
    ψ_g = 1.0/ω * sqrt((V_tI+R*I_I+Xl*I_R)^2 + (V_tR+R*I_R-Xl*I_I)^2)

    # Compute saturation and saturated reactances
    Se = saturation_function(machine, ψ_g)

    Xd_ppsat = Xl + (Xd_pp-Xl)/(1+Se)
    Xq_ppsat = Xl + (Xq_pp-Xl)/(1+Se)


    #Currents
    I_d =
        (1.0 / (R^2 + ω^2 * Xd_ppsat * Xq_ppsat)) *
        (-R * (V_dq[1] + ω*ψq_pp) - ω*Xq_ppsat * (V_dq[2] -ω*ψd_pp))
   
    I_q =
        (1.0 / (R^2 + ω^2 * Xd_ppsat * Xq_ppsat)) * 
        (-R * (V_dq[2] - ω*ψd_pp) + ω*Xd_ppsat * (V_dq[1] +ω*ψq_pp))
    
    Xad_Ifd = (1.0+Se)/(1.0-Kw*I_d) * (eq_p + (Xd - Xd_p) * (I_d/(1.0+Se) + γ_d2 * (eq_p - ψd_p - (Xd_p - Xl)* I_d/(1.0+Se))))

    Xaq_I1q = (-ed_p + (Xq - Xq_p) * (I_q/(1.0+Se) - γ_q2 * (ed_p - ψq_p + (Xq_p - Xl)* I_q/(1.0+Se))))

    # Electric torque
    ψd = ψd_pp - Xd_ppsat*I_d
    ψq = ψq_pp - Xq_ppsat*I_q 

    τ_e = ψd*I_q - ψq*I_d 
    
    #Compute ODEs
    output_ode[local_ix[1]] = (1.0 / Td0_p) * (Vf - Xad_Ifd)                        # eq_p
    output_ode[local_ix[2]] = (1.0+Se) / Tq0_p * (Xaq_I1q)                            # ed_p
    output_ode[local_ix[3]] = (1.0+Se) / Td0_pp * (-ψd_p + eq_p - (Xd_p - Xl) * I_d/(1.0+Se))   # ψd_p
    output_ode[local_ix[4]] = (1.0+Se) / Tq0_pp * (-ψq_p + ed_p + (Xq_p - Xl) * I_q/(1.0+Se))   # ψq_p

    #Update inner_vars
    inner_vars[τe_var] = τ_e
    inner_vars[Xad_Ifd_var] = Xad_Ifd

    #Compute current from the generator to the grid
    I_RI = (basepower / Sbase) * dq_ri(δ) * [I_d; I_q]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end

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
inner_vars,
current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.FullMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.FullMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψd = internal_states[1]
    ψq = internal_states[2]
    ψf = internal_states[3]
    ψ1d = internal_states[4]
    ψ1q = internal_states[5]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.FullMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]
    #ψd = inner_vars[ψd_var]
    #ψq = inner_vars[ψq_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    R_f = PSY.get_R_f(machine)
    R_1d = PSY.get_R_1d(machine)
    R_1q = PSY.get_R_1q(machine)
    inv_d_fluxlink = PSY.get_inv_d_fluxlink(machine)
    inv_q_fluxlink = PSY.get_inv_q_fluxlink(machine)
    basepower = PSY.get_base_power(dynamic_device)

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
    inner_vars[τe_var] = τ_e
    inner_vars[ψd_var] = ψd
    inner_vars[ψq_var] = ψq

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
inner_vars,
current_r::AbstractArray{<:ACCEPTED_REAL_TYPES},
    current_i::AbstractArray{<:ACCEPTED_REAL_TYPES},
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.SimpleFullMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(dynamic_device, PSY.SimpleFullMachine)

    #Define internal states for component
    internal_states = @view device_states[local_ix]
    ψf = internal_states[1]
    ψ1d = internal_states[2]
    ψ1q = internal_states[3]

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(dynamic_device, PSY.SimpleFullMachine)
    δ = device_states[external_ix[1]]
    ω = device_states[external_ix[2]]

    #Obtain inner variables for component
    V_tR = inner_vars[VR_gen_var]
    V_tI = inner_vars[VI_gen_var]
    Vf = inner_vars[Vf_var]
    ψd = inner_vars[ψd_var]
    ψq = inner_vars[ψq_var]

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    inner_vars[τe_var] = τ_e
    inner_vars[ψd_var] = ψd
    inner_vars[ψq_var] = ψq

    #Compute current from the generator to the grid
    I_RI = (BaseMVA / Sbase) * dq_ri(δ) * [i_df1d[1]; i_q1q[1]]

    #Update current
    current_r[1] += I_RI[1]
    current_i[1] += I_RI[2]

    return
end
=#
