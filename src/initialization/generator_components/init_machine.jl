"""
Initialitation of model of 0-state synchronous (classic model) machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(device_states,
    device::PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen) / PSY.get_basepower(static_gen)
    Q0 = PSY.get_reactivepower(static_gen) / PSY.get_basepower(static_gen)
    Vm = PSY.get_voltage(PSY.get_bus(static_gen))
    θ = PSY.get_angle(PSY.get_bus(static_gen))
    S0 = P0 + Q0*1im
    V_R = Vm*cos(θ)
    V_I = Vm*sin(θ)
    V = V_R + V_I*1im
    I = conj(S0/V)

    #Machine Data
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    #Assumption of Classical Machine: Xq = Xd_p
    Xd_p = PSY.get_Xd_p(machine)

    δ0 = angle(V + (R + Xd_p*1im)*I)
    ω0 = 1.0
    τm0 = real(V*conj(I))
    #To solve: δ, τm, Vf0
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]

        V_dq = LITS.ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (Vf0 - V_dq[2]) - R * V_dq[1])  #15.36
        i_q = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_dq[1] + R * (Vf0 - V_dq[2])) #15.36
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q 
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1]*i_d + V_dq[2]*i_q) #Output Power
        out[3] = Q0 - (V_dq[2]*i_d - V_dq[1]*i_q) #Output Reactive Power
    end
    x0 = [δ0, τm0, 1.0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        get_inner_vars(device)[VR_gen_var] = V_R
        get_inner_vars(device)[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        get_inner_vars(device)[τe_var] = sol_x0[2] 
        get_inner_vars(device)[τm_var] = sol_x0[2]
        #Not necessary to update Vf for AVR in Base Machine. Update eq_p:
        PSY.set_eq_p!(machine, sol_x0[3])
        get_inner_vars(device)[Vf_var] = sol_x0[3]
    end
end

"""
Initialitation of model of 2-state (One d- and One q-) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(device_states,
    device::PSY.DynamicGenerator{PSY.OneDOneQMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen) / PSY.get_basepower(static_gen)
    Q0 = PSY.get_reactivepower(static_gen) / PSY.get_basepower(static_gen)
    Vm = PSY.get_voltage(PSY.get_bus(static_gen))
    θ = PSY.get_angle(PSY.get_bus(static_gen))
    S0 = P0 + Q0*1im
    V_R = Vm*cos(θ)
    V_I = Vm*sin(θ)
    V = V_R + V_I*1im
    I = conj(S0/V)

    #Machine Data
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    
    #States of OneDOneQMachine are [1] eq_p and [2] ed_p
    δ0 = angle(V + (R + Xq*1im)*I)
    ω0 = 1.0
    τm0 = real(V*conj(I))
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        eq_p = x[4]
        ed_p = x[5]

        V_dq = LITS.ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / (R^2 + Xd_p * Xq_p)) * (Xq_p * (eq_p - V_dq[2]) + R * (ed_p - V_dq[1]))  #15.32
        i_q = (1.0 / (R^2 + Xd_p * Xq_p)) * (-Xd_p * (ed_p - V_dq[1]) + R * (eq_p - V_dq[2]))  #15.32
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q 
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1]*i_d + V_dq[2]*i_q) #Output Power
        out[3] = Q0 - (V_dq[2]*i_d - V_dq[1]*i_q) #Output Reactive Power
        out[4] = -eq_p - (Xd - Xd_p) * i_d + Vf0 #∂(eq_p)/∂t
        out[5] = -ed_p + (Xq - Xq_p) * i_q #∂(ed_p)/∂t
    end
    V_dq0 = LITS.ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        get_inner_vars(device)[VR_gen_var] = V_R
        get_inner_vars(device)[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        get_inner_vars(device)[τe_var] = sol_x0[2] 
        get_inner_vars(device)[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine. 
        get_inner_vars(device)[Vf_var] = sol_x0[3]
        #Update eq_p and ed_p for Machine
        machine_ix = get_local_state_ix(device, PSY.OneDOneQMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
    end
end

"""
Initialitation of model of 6-state (Marconato) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(device_states,
    device::PSY.DynamicGenerator{PSY.MarconatoMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen) / PSY.get_basepower(static_gen)
    Q0 = PSY.get_reactivepower(static_gen) / PSY.get_basepower(static_gen)
    Vm = PSY.get_voltage(PSY.get_bus(static_gen))
    θ = PSY.get_angle(PSY.get_bus(static_gen))
    S0 = P0 + Q0*1im
    V_R = Vm*cos(θ)
    V_I = Vm*sin(θ)
    V = V_R + V_I*1im
    I = conj(S0/V)

    #Machine Data
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    T_AA = PSY.get_T_AA(machine)
    Td0_p = PSY.get_Td0_p(machine)
    γd = PSY.get_γd(machine)
    γq = PSY.get_γq(machine)
    
    #States of MarconatoMachine are [1] ψq, [2] ψd, [3] eq_p, [4] ed_p, [5] eq_pp and [6] ed_pp
    δ0 = angle(V + (R + Xq*1im)*I)
    ω0 = 1.0
    τm0 = real(V*conj(I))
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        ψq = x[4]
        ψd = x[5]
        eq_p = x[6]
        ed_p = x[7]
        eq_pp = x[8]
        ed_pp = x[9]

        V_dq = LITS.ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / Xd_pp) * (eq_pp - ψd)      #15.18
        i_q = (1.0 / Xq_pp) * (-ed_pp - ψq)     #15.18
        τ_e = ψd * i_q - ψq * i_d               #15.6
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1]*i_d + V_dq[2]*i_q) #Output Power
        out[3] = Q0 - (V_dq[2]*i_d - V_dq[1]*i_q) #Output Reactive Power
        out[4] = R * i_q - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[5] = R * i_d + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[6] = -eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf0    #15.16 eq_p
        out[7] = -ed_p + (Xq - Xq_p - γq) * i_q                                 #15.16 ed_p
        out[8] = -eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf0       #15.16 eq_pp
        out[9] = -ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q #15.16 ed_pp
    end
    
    V_dq0 = LITS.ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[1], V_dq0[2], V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        get_inner_vars(device)[VR_gen_var] = V_R
        get_inner_vars(device)[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        get_inner_vars(device)[τe_var] = sol_x0[2] 
        get_inner_vars(device)[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine. 
        get_inner_vars(device)[Vf_var] = sol_x0[3]
        #Update states for Machine
        machine_ix = get_local_state_ix(device, PSY.MarconatoMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #ψq
        machine_states[2] = sol_x0[5] #ψd
        machine_states[3] = sol_x0[6] #eq_p
        machine_states[4] = sol_x0[7] #ed_p
        machine_states[5] = sol_x0[8] #eq_pp
        machine_states[6] = sol_x0[9] #ed_pp
        #Update fluxes inner vars
        get_inner_vars(device)[ψq_var] = sol_x0[4]
        get_inner_vars(device)[ψd_var] = sol_x0[5]
    end
end



"""
Initialitation of model of 4-state (Simple Marconato) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(device_states,
    device::PSY.DynamicGenerator{PSY.SimpleMarconatoMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen) / PSY.get_basepower(static_gen)
    Q0 = PSY.get_reactivepower(static_gen) / PSY.get_basepower(static_gen)
    Vm = PSY.get_voltage(PSY.get_bus(static_gen))
    θ = PSY.get_angle(PSY.get_bus(static_gen))
    S0 = P0 + Q0*1im
    V_R = Vm*cos(θ)
    V_I = Vm*sin(θ)
    V = V_R + V_I*1im
    I = conj(S0/V)

    #Machine Data
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    T_AA = PSY.get_T_AA(machine)
    Td0_p = PSY.get_Td0_p(machine)
    γd = PSY.get_γd(machine)
    γq = PSY.get_γq(machine)
    
    #States of SimpleMarconatoMachine are [1] eq_p, [2] ed_p, [3] eq_pp and [4] ed_pp
    δ0 = angle(V + (R + Xq*1im)*I)
    ω0 = 1.0
    τm0 = real(V*conj(I))
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        eq_p = x[4]
        ed_p = x[5]
        eq_pp = x[6]
        ed_pp = x[7]

        V_dq = LITS.ri_dq(δ) * [V_R; V_I]
        i_d =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) * (Xq_pp * (eq_pp - V_dq[2]) + R * (ed_pp - V_dq[1]))      #15.25
        i_q =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) * (-Xd_pp * (ed_pp - V_dq[1]) + R * (eq_pp - V_dq[2]))      #15.25 
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q 
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1]*i_d + V_dq[2]*i_q) #Output Power
        out[3] = Q0 - (V_dq[2]*i_d - V_dq[1]*i_q) #Output Reactive Power
        out[4] = -eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf0             #15.16 eq_p
        out[5] = -ed_p + (Xq - Xq_p - γq) * i_q                                         #15.16 ed_p
        out[6] = -eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf0        #15.16 eq_pp
        out[7] = -ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q                              #15.16 ed_pp
    end
    V_dq0 = LITS.ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Synch. Machine failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        get_inner_vars(device)[VR_gen_var] = V_R
        get_inner_vars(device)[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        get_inner_vars(device)[τe_var] = sol_x0[2] 
        get_inner_vars(device)[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine. 
        get_inner_vars(device)[Vf_var] = sol_x0[3]
        #Update eq_p and ed_p for Machine
        machine_ix = get_local_state_ix(device, PSY.SimpleMarconatoMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
        machine_states[3] = sol_x0[6] #eq_pp
        machine_states[4] = sol_x0[7] #ed_pp
    end
end