function initialize_avr!(
    device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = get_inner_vars(device)[Vf_var]
    PSY.set_V_ref!(PSY.get_avr(device), Vf)
    PSY.set_Vf!(PSY.get_avr(device), Vf)
    #Update Control Refs
    device.ext[CONTROL_REFS][V_ref_index] = Vf
end

function initialize_avr!(
    device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf0 = get_inner_vars(device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)
    #Solve V_ref
    V_ref = Vm
    #Set Vf state equals to Vf0
    avr_ix = get_local_state_ix(device, PSY.AVRSimple)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vf0 #δ
    #Set V_ref
    PSY.set_V_ref!(PSY.get_avr(device), V_ref)
    #Update Control Refs
    device.ext[CONTROL_REFS][V_ref_index] = V_ref
end

function initialize_avr!(
    device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = get_inner_vars(device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(device)
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
        device.ext[CONTROL_REFS][V_ref_index] = sol_x0[1]
        #Update AVR states
        avr_ix = get_local_state_ix(device, PSY.AVRTypeI)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #δ
        avr_states[2] = sol_x0[2] #ω
        avr_states[3] = sol_x0[3]
        avr_states[4] = Vm
    end
end

function initialize_avr!(
    device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = get_inner_vars(device)[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(device)
    K0 = PSY.get_K0(avr)
    T1 = PSY.get_T1(avr)
    T2 = PSY.get_T2(avr)
    T3 = PSY.get_T3(avr)
    T4 = PSY.get_T4(avr)
    Tr = PSY.get_Tr(avr)
    Vr_max = PSY.get_Vr_max(avr)
    Vr_min = PSY.get_Vr_min(avr)
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
        device.ext[CONTROL_REFS][V_ref_index] = sol_x0[1]
        #Update AVR states
        avr_ix = get_local_state_ix(device, PSY.AVRTypeII)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #Vf
        avr_states[2] = sol_x0[2] #Vr1
        avr_states[3] = sol_x0[3] #Vr2
        avr_states[4] = Vm #Vm
    end
end
