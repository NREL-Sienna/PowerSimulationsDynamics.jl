function initialize_avr!(device_states,
    device::PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = get_inner_vars(device)[Vf_var]
    PSY.set_V_ref!(PSY.get_avr(device), Vf)
    #Update Control Refs
    device.ext[CONTROL_REFS][V_ref_index] = Vf
end

function initialize_avr!(device_states,
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