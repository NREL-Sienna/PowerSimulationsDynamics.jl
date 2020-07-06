function initialize_machine!(device_states,
    device::PSY.DynamicGenerator{PSY.BaseMachine, S, A, TG, P},
) where {S <: PSY.Shaft, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #PowerFlow Data
    static_gen = PSY.get_static_injector(device)
    P0 = PSY.get_activepower(static_gen)
    Q0 = PSY.get_reactivepower(static_gen)
    Vm = PSY.get_voltage(PSY.get_bus(static_gen))
    θ = PSY.get_angle(PSY.get_bus(static_gen))
    S0 = P0 + Q0*1im
    V_R = Vm*cos(θ)
    V_I = Vm*sin(θ)
    V = V_R + V_I*1im
    I = conj(S0/V)
    I_R = real(I)
    I_I = imag(I)

    #Machine Data
    machine = PSY.get_machine(device)
    R = PSY.get_R(machine)
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