"""
Initialitation of model of 0-state synchronous (classic model) machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    #Assumption of Classical Machine: Xq = Xd_p
    Xd_p = PSY.get_Xd_p(machine)

    δ0 = angle(V + (R + Xd_p * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE) τm0, P0
    #To solve: δ, τm, Vf0
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (Vf0 - V_dq[2]) - R * V_dq[1])  #15.36
        i_q = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_dq[1] + R * (Vf0 - V_dq[2])) #15.36
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
    end
    x0 = [δ0, τm0, 1.0]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Not necessary to update Vf for AVR in Base Machine. Update eq_p:
        PSY.set_eq_p!(machine, sol_x0[3])
        inner_vars[Vf_var] = sol_x0[3]
    end
    return
end

"""
Initialitation of model of 2-state (One d- and One q-) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.OneDOneQMachine, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)

    #States of OneDOneQMachine are [1] eq_p and [2] ed_p
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE) τm0, P0
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        eq_p = x[4]
        ed_p = x[5]

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / (R^2 + Xd_p * Xq_p)) * (Xq_p * (eq_p - V_dq[2]) + R * (ed_p - V_dq[1]))  #15.32
        i_q =
            (1.0 / (R^2 + Xd_p * Xq_p)) * (-Xd_p * (ed_p - V_dq[1]) + R * (eq_p - V_dq[2]))  #15.32
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = -eq_p - (Xd - Xd_p) * i_d + Vf0 #∂(eq_p)/∂t
        out[5] = -ed_p + (Xq - Xq_p) * i_q #∂(ed_p)/∂t
    end
    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update eq_p and ed_p for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.OneDOneQMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
    end
    return
end

"""
Initialitation of model of 6-state (SauerPai) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.SauerPaiMachine, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    Xl = PSY.get_Xl(machine)
    Td0_p = PSY.get_Td0_p(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    γ_q2 = PSY.get_γ_q2(machine)

    #States of SauerPaiMachine are [1] ψq, [2] ψd, [3] eq_p, [4] ed_p, [5] ψd_pp and [6] ψq_pp
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE)
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        ψq = x[4]
        ψd = x[5]
        eq_p = x[6]
        ed_p = x[7]
        ψd_pp = x[8]
        ψq_pp = x[9]

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / Xd_pp) * (γ_d1 * eq_p - ψd + (1 - γ_d1) * ψd_pp)      #15.15
        i_q = (1.0 / Xq_pp) * (-γ_q1 * ed_p - ψq + (1 - γ_q1) * ψq_pp)     #15.15
        τ_e = ψd * i_q - ψq * i_d               #15.6
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = R * i_q - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[5] = R * i_d + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[6] = -eq_p - (Xd - Xd_p) * (i_d - γ_d2 * ψd_pp - (1 -γ_d1) * i_d + γ_d2 * eq_p) + Vf0    #15.13
        out[7] = -ed_p + (Xq - Xq_p) * (i_q - γ_q2 * ψq_pp - (1 -γ_q1) * i_q - γ_d2 * ed_p)          #15.13
        out[8] = -ψd_pp + eq_p - (Xd_p - Xl) * i_d      #15.13
        out[9] = -ψq_pp - ed_p - (Xq_p - Xl) * i_q #15.13
    end

    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[1], V_dq0[2], V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.SauerPaiMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #ψq
        machine_states[2] = sol_x0[5] #ψd
        machine_states[3] = sol_x0[6] #eq_p
        machine_states[4] = sol_x0[7] #ed_p
        machine_states[5] = sol_x0[8] #eq_pp
        machine_states[6] = sol_x0[9] #ed_pp
        #Update fluxes inner vars
        inner_vars[ψq_var] = sol_x0[4]
        inner_vars[ψd_var] = sol_x0[5]
    end
    return
end

"""
Initialitation of model of 6-state (Marconato) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.MarconatoMachine, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
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
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE)
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

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / Xd_pp) * (eq_pp - ψd)      #15.18
        i_q = (1.0 / Xq_pp) * (-ed_pp - ψq)     #15.18
        τ_e = ψd * i_q - ψq * i_d               #15.6
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = R * i_q - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[5] = R * i_d + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[6] = -eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf0    #15.16 eq_p
        out[7] = -ed_p + (Xq - Xq_p - γq) * i_q                                 #15.16 ed_p
        out[8] = -eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf0       #15.16 eq_pp
        out[9] = -ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q #15.16 ed_pp
    end

    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[1], V_dq0[2], V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.MarconatoMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #ψq
        machine_states[2] = sol_x0[5] #ψd
        machine_states[3] = sol_x0[6] #eq_p
        machine_states[4] = sol_x0[7] #ed_p
        machine_states[5] = sol_x0[8] #eq_pp
        machine_states[6] = sol_x0[9] #ed_pp
        #Update fluxes inner vars
        inner_vars[ψq_var] = sol_x0[4]
        inner_vars[ψd_var] = sol_x0[5]
    end
    return
end

"""
Initialitation of model of 4-state (Simple Marconato) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SimpleMarconatoMachine, S, A, TG, P},
    },
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
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
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE)
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        eq_p = x[4]
        ed_p = x[5]
        eq_pp = x[6]
        ed_pp = x[7]

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (Xq_pp * (eq_pp - V_dq[2]) + R * (ed_pp - V_dq[1]))      #15.25
        i_q =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (-Xd_pp * (ed_pp - V_dq[1]) + R * (eq_pp - V_dq[2]))      #15.25
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = -eq_p - (Xd - Xd_p - γd) * i_d + (1 - (T_AA / Td0_p)) * Vf0             #15.16 eq_p
        out[5] = -ed_p + (Xq - Xq_p - γq) * i_q                                         #15.16 ed_p
        out[6] = -eq_pp + eq_p - (Xd_p - Xd_pp + γd) * i_d + (T_AA / Td0_p) * Vf0        #15.16 eq_pp
        out[7] = -ed_pp + ed_p + (Xq_p - Xq_pp + γq) * i_q                              #15.16 ed_pp
    end
    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update eq_p and ed_p for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.SimpleMarconatoMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
        machine_states[3] = sol_x0[6] #eq_pp
        machine_states[4] = sol_x0[7] #ed_pp
    end
    return
end

"""
Initialitation of model of 6-state (Anderson-Fouad) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.AndersonFouadMachine, S, A, TG, P},
    },
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)

    #States of Anderson-Fouad are [1] ψq, [2] ψd, [3] eq_p, [4] ed_p, [5] eq_pp and [6] ed_pp
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))

    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE)
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

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d = (1.0 / Xd_pp) * (eq_pp - ψd)      #15.18
        i_q = (1.0 / Xq_pp) * (-ed_pp - ψq)     #15.18
        τ_e = ψd * i_q - ψq * i_d               #15.6
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = R * i_q - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[5] = R * i_d + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[6] = -eq_p + (Xd - Xd_p) * i_d + Vf0           #15.19 eq_p
        out[7] = -ed_p + (Xq - Xq_p) * i_q                 #15.19 ed_p
        out[8] = -eq_pp + eq_p - (Xd_p - Xd_pp) * i_d      #15.19 eq_pp
        out[9] = -ed_pp + ed_p + (Xq_p - Xq_pp) * i_q      #15.19 ed_pp
    end

    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[1], V_dq0[2], V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.AndersonFouadMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #ψq
        machine_states[2] = sol_x0[5] #ψd
        machine_states[3] = sol_x0[6] #eq_p
        machine_states[4] = sol_x0[7] #ed_p
        machine_states[5] = sol_x0[8] #eq_pp
        machine_states[6] = sol_x0[9] #ed_pp
        #Update fluxes inner vars
        inner_vars[ψq_var] = sol_x0[4]
        inner_vars[ψd_var] = sol_x0[5]
    end
    return
end

"""
Initialitation of model of 4-state (Simple Anderson-Fouad) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{PSY.SimpleAFMachine, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xq = PSY.get_Xq(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)

    #States of SimpleAndersonFouad are [1] eq_p, [2] ed_p, [3] eq_pp and [4] ed_pp
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))
    @assert isapprox(τm0, P0; atol = STRICT_NLSOLVE_F_TOLERANCE)
    #To solve: δ, τm, Vf0, eq_p, ed_p
    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        eq_p = x[4]
        ed_p = x[5]
        eq_pp = x[6]
        ed_pp = x[7]

        V_dq = ri_dq(δ) * [V_R; V_I]
        i_d =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (Xq_pp * (eq_pp - V_dq[2]) + R * (ed_pp - V_dq[1]))      #15.25
        i_q =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (-Xd_pp * (ed_pp - V_dq[1]) + R * (eq_pp - V_dq[2]))      #15.25
        Pe = (V_dq[1] + R * i_d) * i_d + (V_dq[2] + R * i_q) * i_q
        out[1] = τm - Pe #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_d + V_dq[2] * i_q) #Output Power
        out[3] = Q0 - (V_dq[2] * i_d - V_dq[1] * i_q) #Output Reactive Power
        out[4] = -eq_p + (Xd - Xd_p) * i_d + Vf0           #15.19 eq_p
        out[5] = -ed_p + (Xq - Xq_p) * i_q                 #15.19 ed_p
        out[6] = -eq_pp + eq_p - (Xd_p - Xd_pp) * i_d      #15.19 eq_pp
        out[7] = -ed_pp + ed_p + (Xq_p - Xq_pp) * i_q      #15.19 ed_pp
    end
    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], V_dq0[1], V_dq0[2], V_dq0[1]]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update eq_p and ed_p for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.SimpleAFMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
        machine_states[3] = sol_x0[6] #eq_pp
        machine_states[4] = sol_x0[7] #ed_pp
    end
    return
end

function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, A, TG, P}},
    inner_vars::AbstractVector,
) where {
    M <: Union{PSY.RoundRotorQuadratic, PSY.RoundRotorExponential},
    S <: PSY.Shaft,
    A <: PSY.AVR,
    TG <: PSY.TurbineGov,
    P <: PSY.PSS,
}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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
    Se = PSY.get_Se(machine)
    # Initialization doesn't consider saturation
    #Sat_A, Sat_B = PSY.get_saturation_coeffs(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)
    γ_q2 = PSY.get_γ_q2(machine)
    γ_qd = PSY.get_γ_qd(machine)

    ## Initialization ##
    ## Fluxes
    ψ_pp0 = V + (R + Xd_pp * 1im) * I
    ψ_pp0_ang = angle(ψ_pp0)
    ψ_pp0_abs = abs(ψ_pp0)
    Se0 = saturation_function(machine, ψ_pp0_abs)
    ## Angles
    _a = ψ_pp0_abs * (Se0 * γ_qd + 1.0)
    _b = (Xq_pp - Xq) * abs(I)
    θ_It = ψ_pp0_ang - angle(I)
    δ0 = ψ_pp0_ang + atan((_b * cos(θ_It)) / (_b * sin(θ_It) - _a))
    #Fluxes
    ψ_pp0_dq = ri_dq(δ0 + pi / 2) * [real(ψ_pp0); imag(ψ_pp0)]
    ψd_pp0 = ψ_pp0_dq[1]
    ψq_pp0 = -ψ_pp0_dq[2]
    ## Currents
    I_d0, I_q0 = ri_dq(δ0) * [real(I); imag(I)]
    ## Voltages and Fluxes
    V_d0, V_q0 = ri_dq(δ0) * [V_R; V_I]
    ψq_pp0 = V_d0 - R * I_d0 - Xq_pp * I_q0
    ψd_pp0 = V_q0 + I_d0 * Xd_pp + I_q0 * R
    # States
    ## External Variables
    τm0 = I_d0 * (V_d0 + I_d0 * R) + I_q0 * (V_q0 + I_q0 * R)
    Vf0 = I_d0 * (Xd - Xd_pp) + ψd_pp0 * (Se0 + 1.0)
    ## States
    eq_p0 = I_d0 * (Xd_p - Xd) - Se0 * ψd_pp0 + Vf0
    ed_p0 = I_q0 * (Xq - Xq_p) - Se0 * γ_qd * ψq_pp0
    ψ_kd0 = -I_d0 * (Xd - Xl) - Se0 * ψd_pp0 + Vf0
    ψ_kq0 = I_q0 * (Xq - Xl) - Se0 * γ_qd * ψq_pp0
    #Auxiliar Variables
    Xad_Ifd0 =
        ed_p0 +
        (Xq - Xq_p) * (γ_q2 * ed_p0 - γ_q2 * ψ_kq0 - γ_q1 * I_q0) +
        Se0 * ψq_pp0 * γ_qd

    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf = x[3]
        eq_p = x[4]
        ed_p = x[5]
        ψ_kd = x[6]
        ψ_kq = x[7]
        Xad_Ifd_aux = x[8]

        V_dq = ri_dq(δ) * [V_R; V_I]
        ψq_pp = γ_q1 * ed_p + ψ_kq * (1 - γ_q1)
        ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψ_kd
        ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)
        I_d =
            (1.0 / (R^2 + Xq_pp * Xd_pp)) *
            (-R * (V_dq[1] - ψq_pp) + Xq_pp * (-V_dq[2] + ψd_pp))
        I_q =
            (1.0 / (R^2 + Xq_pp * Xd_pp)) *
            (Xd_pp * (V_dq[1] - ψq_pp) + R * (-V_dq[2] + ψd_pp))
        Se = saturation_function(machine, ψ_pp)
        Xad_Ifd = eq_p + (Xd - Xd_p) * (γ_d1 * I_d - γ_d2 * ψ_kd + γ_d2 * eq_p) + Se * ψd_pp
        Xaq_I1q =
            ed_p +
            (Xq - Xq_p) * (γ_q2 * ed_p - γ_q2 * ψ_kq - γ_q1 * I_q) +
            Se * ψq_pp * γ_qd
        τ_e = I_d * (V_dq[1] + I_d * R) + I_q * (V_dq[2] + I_q * R)

        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1] * I_d + V_dq[2] * I_q) #Output Power
        out[3] = Q0 - (V_dq[2] * I_d - V_dq[1] * I_q) #Output Reactive Power
        out[4] = (1.0 / Td0_p) * (Vf - Xad_Ifd) #deq_p/dt
        out[5] = (1.0 / Tq0_p) * (-Xaq_I1q) #ded_p/dt
        out[6] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d) #dψ_kd/dt
        out[7] = (1.0 / Tq0_pp) * (-ψ_kq + ed_p + (Xq_p - Xl) * I_q) #deq_pp/dt
        out[8] = Xad_Ifd_aux - Xad_Ifd
    end
    x0 = [δ0, τm0, Vf0, eq_p0, ed_p0, ψ_kd0, ψ_kq0, Xad_Ifd0]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = 1.0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in GENROU Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update Xad_Ifd for AVR in GENROU Machine
        inner_vars[Xad_Ifd_var] = sol_x0[8]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, typeof(machine))
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ed_p
        machine_states[3] = sol_x0[6] #ψ_kd
        machine_states[4] = sol_x0[7] #ψ_kq
    end
    return
end

function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleQuadratic, S, A, TG, P},
    },
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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

    ## Initialization ##
    E = V + (R + Xq * 1im) * I
    δ0 = angle(E)
    ## Currents
    I_d0, I_q0 = ri_dq(δ0) * [real(I); imag(I)]
    ## Voltages and fluxes
    V_d0, V_q0 = ri_dq(δ0) * [V_R; V_I]
    ψ_d0 = V_q0 + R * I_q0
    ψ_q0 = -V_d0 - R * I_d0
    ψd_pp0 = V_q0 + I_d0 * Xd_pp + I_q0 * R
    # States
    ψq_pp0 = -(V_d0 - I_q0 * Xq_pp - I_d0 * R)
    eq_p0 = (1 / (γ_d1 + γ_q1)) * (ψd_pp0 + I_d0 * (Xd_p - Xl) * γ_q1)
    ψ_kd0 = eq_p0 - I_d0 * (Xd_p - Xl)
    ## External Variables
    τm0 = ψ_d0 * I_q0 - ψ_q0 * I_d0
    Se0 = saturation_function(machine, eq_p0)
    Vf0 = eq_p0 + I_d0 * (Xd - Xd_p) + Se0 * eq_p0
    Xad_Ifd0 =
        eq_p0 +
        Se0 * eq_p0 +
        (Xd - Xd_p) * (I_d0 + γ_d2 * (eq_p0 - ψ_kd0 - (Xd_p - Xl) * I_d0))

    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf = x[3]
        eq_p = x[4]
        ψ_kd = x[5]
        ψq_pp = x[6]
        Xad_Ifd_aux = x[7]

        V_d, V_q = ri_dq(δ) * [V_R; V_I]
        #Additional Fluxes
        ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd

        #Currents
        I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d + ψq_pp) + Xd_pp * (ψd_pp - V_q))
        I_q = (1.0 / (R^2 + Xd_pp^2)) * (Xd_pp * (V_d + ψq_pp) + R * (ψd_pp - V_q))
        Se = saturation_function(machine, eq_p)
        Xad_Ifd =
            eq_p +
            Se * eq_p +
            (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))
        τ_e = I_d * (V_d + I_d * R) + I_q * (V_q + I_q * R)

        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_d * I_d + V_q * I_q) #Output Power
        out[3] = Q0 - (V_q * I_d - V_d * I_q) #Output Reactive Power
        out[4] = (1.0 / Td0_p) * (Vf - Xad_Ifd) #deq_p/dt
        out[5] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d) #dψ_kd/dt
        out[6] = (1.0 / Tq0_pp) * (-ψq_pp - (Xq - Xq_pp) * I_q) #dψq_pp/dt
        out[7] = Xad_Ifd_aux - Xad_Ifd
    end
    x0 = [δ0, τm0, Vf0, eq_p0, ψ_kd0, ψq_pp0, Xad_Ifd0]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = 1.0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in GENSAL Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update Xad_Ifd for AVR in GENSAL Machine
        inner_vars[Xad_Ifd_var] = sol_x0[7]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, typeof(machine))
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ψ_kd
        machine_states[3] = sol_x0[6] #ψq_pp
    end
    return
end

function initialize_mach_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{
        PSY.DynamicGenerator{PSY.SalientPoleExponential, S, A, TG, P},
    },
    inner_vars::AbstractVector,
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Get parameters
    machine = PSY.get_machine(dynamic_device)
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

    ## Initialization ##
    ## Fluxes
    ψ_pp0 = V + (R + Xd_pp * 1im) * I
    ψ_pp0_ang = angle(ψ_pp0)
    ψ_pp0_abs = abs(ψ_pp0)
    Se0 = saturation_function(machine, ψ_pp0_abs)
    ## Angles
    _a = ψ_pp0_abs * (Se0 * γ_qd + 1.0)
    _b = (Xq_pp - Xq) * abs(I)
    θ_It = ψ_pp0_ang - angle(I)
    δ0 = ψ_pp0_ang + atan((_b * cos(θ_It)) / (_b * sin(θ_It) - _a))
    #Fluxes
    ψ_pp0_dq = ri_dq(δ0 + pi / 2) * [real(ψ_pp0); imag(ψ_pp0)]
    ψd_pp0 = ψ_pp0_dq[1]
    ψq_pp0 = -ψ_pp0_dq[2]
    ## Currents
    I_d0, I_q0 = ri_dq(δ0) * [real(I); imag(I)]
    ## Voltages and Fluxes
    V_d0, V_q0 = ri_dq(δ0) * [V_R; V_I]
    ψq_pp0 = V_d0 - R * I_d0 - Xq_pp * I_q0
    ψd_pp0 = V_q0 + I_d0 * Xd_pp + I_q0 * R
    ψ_d0 = V_q0 + R * I_q0
    ψ_q0 = -V_d0 - R * I_d0
    # States
    eq_p0 = (1 / (γ_d1 + γ_q1)) * (ψd_pp0 + I_d0 * (Xd_p - Xl) * γ_q1)
    ψ_kd0 = eq_p0 - I_d0 * (Xd_p - Xl)
    ## External Variables
    τm0 = ψ_d0 * I_q0 - ψ_q0 * I_d0
    Vf0 = eq_p0 + I_d0 * (Xd - Xd_p) + ψd_pp0 * Se0
    Xad_Ifd0 =
        eq_p0 +
        Se0 * ψd_pp0 +
        (Xd - Xd_p) * (I_d0 + γ_d2 * (eq_p0 - ψ_kd0 - (Xd_p - Xl) * I_d0))

    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf = x[3]
        eq_p = x[4]
        ψ_kd = x[5]
        ψq_pp = x[6]
        Xad_Ifd_aux = x[7]

        V_d, V_q = ri_dq(δ) * [V_R; V_I]
        #Additional Fluxes
        ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd
        ψ_pp = sqrt(ψd_pp^2 + ψq_pp^2)

        I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d - ψq_pp) + Xq_pp * (-V_q + ψd_pp))
        I_q = (1.0 / (R^2 + Xd_pp^2)) * (Xd_pp * (V_d - ψq_pp) + R * (-V_q + ψd_pp))
        Se = saturation_function(machine, ψ_pp)
        Xad_Ifd =
            eq_p +
            Se * ψd_pp +
            (Xd - Xd_p) * (I_d + γ_d2 * (eq_p - ψ_kd - (Xd_p - Xl) * I_d))
        τ_e = I_d * (V_d + I_d * R) + I_q * (V_q + I_q * R)

        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_d * I_d + V_q * I_q) #Output Power
        out[3] = Q0 - (V_q * I_d - V_d * I_q) #Output Reactive Power
        out[4] = (1.0 / Td0_p) * (Vf - Xad_Ifd) #deq_p/dt
        out[5] = (1.0 / Td0_pp) * (-ψ_kd + eq_p - (Xd_p - Xl) * I_d) #dψ_kd/dt
        out[6] = (1.0 / Tq0_pp) * (-ψq_pp + (Xq - Xq_pp) * I_q - Se * γ_qd * ψq_pp) #dψq_pp/dt
        out[7] = Xad_Ifd_aux - Xad_Ifd
    end
    x0 = [δ0, τm0, Vf0, eq_p0, ψ_kd0, ψq_pp0, Xad_Ifd0]
    sol = NLsolve.nlsolve(f!, x0, ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        inner_vars[VR_gen_var] = V_R
        inner_vars[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = 1.0 #ω
        #Update Mechanical and Electrical Torque on Generator
        inner_vars[τe_var] = sol_x0[2]
        inner_vars[τm_var] = sol_x0[2]
        #Update Vf for AVR in GENSAL Machine.
        inner_vars[Vf_var] = sol_x0[3]
        #Update Xad_Ifd for AVR in GENSAL Machine
        inner_vars[Xad_Ifd_var] = sol_x0[7]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, typeof(machine))
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #eq_p
        machine_states[2] = sol_x0[5] #ψ_kd
        machine_states[3] = sol_x0[6] #ψq_pp
    end
    return
end

#=
"""
Initialitation of model of 5-state (Kundur) synchronous machine in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""
function initialize_mach_shaft!(
    device_states,
static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{PSY.FullMachine, S, A, TG, P}},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data

    P0 = PSY.get_active_power(static)
    Q0 = PSY.get_reactive_power(static)
    Vm = PSY.get_magnitude(PSY.get_bus(static))
    θ = PSY.get_angle(PSY.get_bus(static))
    S0 = P0 + Q0 * 1im
    V_R = Vm * cos(θ)
    V_I = Vm * sin(θ)
    V = V_R + V_I * 1im
    I = conj(S0 / V)

    #Machine Data
    machine = PSY.get_machine(dynamic_device)
    R = PSY.get_R(machine)
    R_f = PSY.get_R_f(machine)
    R_1d = PSY.get_R_1d(machine)
    R_1q = PSY.get_R_1q(machine)
    Xq = PSY.get_L_q(machine)'
    Xd = PSY.get_L_d(machine)
    Xad = PSY.get_L_ad(machine)
    inv_d_fluxlink = PSY.get_inv_d_fluxlink(machine)
    inv_q_fluxlink = PSY.get_inv_q_fluxlink(machine)

    #States of FullMachine are [1] ψd, [2] ψq, [3] ψf , [4] ψ1d, [5] ψ1q
    δ0 = angle(V + (R + Xq * 1im) * I)
    ω0 = 1.0
    τm0 = real(V * conj(I))

    function f!(out, x)
        δ = x[1]
        τm = x[2]
        Vf0 = x[3]
        ψd = x[4]
        ψq = x[5]
        ψf = x[6]
        ψ1d = x[7]
        ψ1q = x[8]

        V_dq = ri_dq(δ) * [V_R; V_I]
        #Obtain electric currents] using Flux Linkage equations
        i_df1d = inv_d_fluxlink * [ψd; ψf; ψ1d]  #11.18 in Machowski (3.127, 3.130 and 3.131 in Kundur)
        i_q1q = inv_q_fluxlink * [ψq; ψ1q]        #11.19 in Machowski (3.128 and 3.132 in Kundur)
        τ_e = ψd * i_q1q[1] - ψq * i_df1d[1]            #15.6 in Milano or 3.117 in Kundur
        out[1] = τm - τ_e #Mechanical Torque
        out[2] = P0 - (V_dq[1] * i_df1d[1] + V_dq[2] * i_q1q[1]) #Output Power
        out[3] = Q0 - (V_dq[2] * i_df1d[1] - V_dq[1] * i_q1q[1]) #Output Reactive Power
        out[4] = R * i_df1d[1] + ω0 * ψq + V_dq[1]                                    #15.9 ψd
        out[5] = R * i_q1q[1] - ω0 * ψd + V_dq[2]                                    #15.9 ψq
        out[6] = -R_f * i_df1d[2] + Vf0          #15.19 eq_p
        out[7] = -R_1d * i_df1d[3]                #15.19 ed_p
        out[8] = -R_1q * i_q1q[2]        #15.19 eq_pp
    end

    V_dq0 = ri_dq(δ0) * [V_R; V_I]
    x0 = [δ0, τm0, 1.0, V_dq0[2], -V_dq0[1], 1.0, 0.0, 0.0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Machine $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update terminal voltages
        get_inner_vars(dynamic_device)[VR_gen_var] = V_R
        get_inner_vars(dynamic_device)[VI_gen_var] = V_I
        #Update δ and ω of Shaft. Works for every Shaft.
        shaft_ix = get_local_state_ix(dynamic_device, S)
        shaft_states = @view device_states[shaft_ix]
        shaft_states[1] = sol_x0[1] #δ
        shaft_states[2] = ω0 #ω
        #Update Mechanical and Electrical Torque on Generator
        get_inner_vars(dynamic_device)[τe_var] = sol_x0[2]
        get_inner_vars(dynamic_device)[τm_var] = sol_x0[2]
        #Update Vf for AVR in OneDOneQ Machine.
        get_inner_vars(dynamic_device)[Vf_var] = sol_x0[3]
        #Update states for Machine
        machine_ix = get_local_state_ix(dynamic_device, PSY.FullMachine)
        machine_states = @view device_states[machine_ix]
        machine_states[1] = sol_x0[4] #ψd
        machine_states[2] = sol_x0[5] #ψq
        machine_states[3] = sol_x0[6] #ψf
        machine_states[4] = sol_x0[7] #ψ1d
        machine_states[5] = sol_x0[8] #ψ1q
        #Update fluxes inner vars
        get_inner_vars(dynamic_device)[ψd_var] = sol_x0[4]
        get_inner_vars(dynamic_device)[ψq_var] = sol_x0[5]
    end
end
=#
