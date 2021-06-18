function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = get_inner_vars(dynamic_device)[Vf_var]
    PSY.set_V_ref!(PSY.get_avr(dynamic_device), Vf)
    PSY.set_Vf!(PSY.get_avr(dynamic_device), Vf)
    #Update Control Refs
    PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = Vf
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf0 = get_inner_vars(dynamic_device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(
        get_inner_vars(dynamic_device)[VR_gen_var]^2 +
        get_inner_vars(dynamic_device)[VI_gen_var]^2,
    )
    #Solve V_ref
    V_ref = Vm
    #Set Vf state equals to Vf0
    avr_ix = get_local_state_ix(dynamic_device, PSY.AVRSimple)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vf0 #δ
    #Set V_ref
    PSY.set_V_ref!(PSY.get_avr(dynamic_device), V_ref)
    #Update Control Refs
    PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = V_ref
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = get_inner_vars(dynamic_device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(
        get_inner_vars(dynamic_device)[VR_gen_var]^2 +
        get_inner_vars(dynamic_device)[VI_gen_var]^2,
    )

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ka = PSY.get_Ka(avr)
    Ke = PSY.get_Ke(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)
    #Obtain saturated Vf
    Se_Vf = Ae * exp(Be * abs(Vf0))

    #States of AVRTypeI are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr1, Vr2
    function f!(out, x)
        V_ref = x[1]
        Vr1 = x[2]
        Vr2 = x[3]

        out[1] = Vf0 * (Ke + Se_Vf) - Vr1 #16.12c
        out[2] = Ka * (V_ref - Vm - Vr2 - (Kf / Tf) * Vf0) - Vr1 #16.12a
        out[3] = (Kf / Tf) * Vf0 + Vr2 #16.12b
    end
    x0 = [1.0, Vf0, Vf0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in AVR failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = sol_x0[1]
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeI)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #δ
        avr_states[2] = sol_x0[2] #ω
        avr_states[3] = sol_x0[3]
        avr_states[4] = Vm
    end
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = get_inner_vars(dynamic_device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(
        get_inner_vars(dynamic_device)[VR_gen_var]^2 +
        get_inner_vars(dynamic_device)[VI_gen_var]^2,
    )

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    K0 = PSY.get_K0(avr)
    T1 = PSY.get_T1(avr)
    T2 = PSY.get_T2(avr)
    T3 = PSY.get_T3(avr)
    T4 = PSY.get_T4(avr)
    Tr = PSY.get_Tr(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)
    #Obtain saturated Vf
    Se_Vf = Ae * (exp(Be * abs(Vf0)))

    #States of AVRTypeII are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr1, Vr2
    function f!(out, x)
        V_ref = x[1]
        Vr1 = x[2]
        Vr2 = x[3]

        Vr = K0 * Vr2 + (T4 / T3) * (Vr1 + K0 * (T2 / T1) * (V_ref - Vm)) #16.21

        out[1] = (Vf0 * (1.0 + Se_Vf) - Vr) #16.18
        out[2] = K0 * (1.0 - (T2 / T1)) * (V_ref - Vm) - Vr1 #16.14
        out[3] = (1.0 - (T4 / T3)) * (Vr1 + K0 * (T2 / T1) * (V_ref - Vm)) - K0 * Vr2  #16.20
    end
    x0 = [1.0, Vf0, Vf0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in AVR failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = sol_x0[1]
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeII)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #Vf
        avr_states[2] = sol_x0[2] #Vr1
        avr_states[3] = sol_x0[3] #Vr2
        avr_states[4] = Vm #Vm
    end
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::PSY.DynamicGenerator{M, S, PSY.ESAC1A, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = get_inner_vars(dynamic_device)[Vf_var]
    #Obtain I_fd obtained from Machine:
    Xad_Ifd0 = get_inner_vars(dynamic_device)[Xad_Ifd_var]
    #Obtain measured terminal voltage
    Vm0 = sqrt(
        get_inner_vars(dynamic_device)[VR_gen_var]^2 +
        get_inner_vars(dynamic_device)[VI_gen_var]^2,
    )

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Tr = PSY.get_Tr(avr)
    Tb = PSY.get_Tb(avr)
    Tc = PSY.get_Tc(avr)
    Ka = PSY.get_Ka(avr)
    Ta = PSY.get_Ta(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    Te = PSY.get_Te(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    Kc = PSY.get_Kc(avr)
    Kd = PSY.get_Kd(avr)
    Ke = PSY.get_Ke(avr)
    E1, E2 = PSY.get_E_sat(avr)
    SE1, SE2 = PSY.get_Se(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    #Obtain saturation
    #Se_Vf = saturation_function(Vm)

    #Solve Ve from rectifier function
    function f_Ve!(out, x)
        V_e0 = x[1]
        I_N0 = Kc * Xad_Ifd0 / V_e0

        out[1] = Vf0 - V_e0 * rectifier_function(I_N0)
    end
    x0 = [1.0]
    sol = NLsolve.nlsolve(f_Ve!, x0)
    if !NLsolve.converged(sol)
        @error("Initialization in V_e0 failed")
    else
        sol_x0 = sol.zero
        V_e0 = sol_x0[1]
        I_N0 = Kc * Xad_Ifd0 / V_e0
    end
    Se0 = saturation_function(avr, V_e0) #To do saturation function
    V_FE0 = Kd * Xad_Ifd0 + Ke * V_e0 + Se0 * V_e0
    V_r20 = V_FE0
    if (V_r20 > Vr_max) || (V_r20 < Vr_min)
        @error("Regulator Voltage V_R = $(V_r20) outside the limits")
    end
    Tc_Tb_ratio = Tb <= eps() ? 0.0 : Tc / Tb
    V_r30 = -(Kf / Tf) * V_FE0
    V_r10 = (V_r20 / Ka) - (Tc_Tb_ratio) * (V_r20 / Ka)
    V_ref0 = Vm0 + V_r20 / Ka

    #States of ESAC1A are Vm, Vr1, Vr2, Ve, Vr3
    #To solve V_ref, Vr1, Vr2, Ve, Vr3
    function f!(out, x)
        V_ref = x[1]
        Vr1 = x[2]
        Vr2 = x[3]
        Ve = x[4]
        Vr3 = x[5]

        #Compute auxiliary variables

        I_N = Kc * Xad_Ifd0 / Ve
        Se = saturation_function(avr, Ve)
        V_FE = Kd * Xad_Ifd0 + Ke * Ve + Ve * Se
        V_F = Vr3 + (Kf / Tf) * V_FE
        V_in = V_ref - Vm0 - V_F
        V_out = Vr1 + (Tc_Tb_ratio) * V_in

        # Time Constants eliminated because don't appear in SteadyState
        out[1] = (V_in * (1 - Tc_Tb_ratio) - Vr1) #dVr1/dt
        out[2] = (Ka * V_out - Vr2) #dVr2/dt
        out[3] = (-(Kf / Tf) * V_FE - Vr3) #dVr3/dt
        out[4] = (Vr2 - V_FE) #dVe/dt
        out[5] = Vf0 - Ve * rectifier_function(I_N)
    end
    x0 = [V_ref0, V_r10, V_r20, V_e0, V_r30]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in AVR failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        PSY.get_ext(dynamic_device)[CONTROL_REFS][V_ref_index] = sol_x0[1]
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, typeof(avr))
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vm0
        avr_states[2] = sol_x0[2] #Vr1
        avr_states[3] = sol_x0[3] #Vr2
        avr_states[4] = sol_x0[4] #Ve
        avr_states[5] = sol_x0[5] #Vr3
    end
end