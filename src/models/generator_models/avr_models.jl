function mdl_avr_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #TODO Change EMF name for Vf in PowerSystems
    #Update Vf voltage on inner vars
    get_inner_vars(device)[Vf_var] = PSY.get_Vf(PSY.get_avr(device))

    return
end

function mdl_avr_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.AVRSimple)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]

    #Define external states for device
    V_th = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)

    #Get Parameters
    Kv = PSY.get_Kv(PSY.get_avr(device))

    #Compute ODEs
    output_ode[local_ix[1]] = Kv * (V_ref - V_th)

    #Update inner_vars
    get_inner_vars(device)[Vf_var] = Vf

    return
end

"""
Model of AVR Type I in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_avr_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.AVRTypeI)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm = internal_states[4]

    #Define external states for device
    V_th = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)
    Vs = get_inner_vars(device)[V_pss_var]

    #Get parameters
    avr = PSY.get_avr(device)
    Ka = PSY.get_Ka(avr)
    Ke = PSY.get_Ke(avr)
    Kf = PSY.get_Kf(avr)
    Ta = PSY.get_Ta(avr)
    Te = PSY.get_Te(avr)
    Tf = PSY.get_Tf(avr)
    Tr = PSY.get_Tr(avr)
    Vr_max = PSY.get_Vr_max(avr)
    Vr_min = PSY.get_Vr_min(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)

    #Compute auxiliary parameters
    Se_Vf = Ae * exp(Be * abs(Vf)) #16.13
    V_ref = V0_ref + Vs

    #Set anti-windup for Vr1. #TODO in callbacks
    #if Vr1 > Vr_max
    #    Vr1 = Vr_max
    #elseif Vr1 < Vr_min
    #    Vr1 = Vr_min
    #end

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = -(1.0 / Te) * (Vf * (Ke + Se_Vf) - Vr1) #16.12c
    output_ode[local_ix[2]] = (1.0 / Ta) * (Ka * (V_ref - Vm - Vr2 - (Kf / Tf) * Vf) - Vr1) #16.12a
    output_ode[local_ix[3]] = -(1.0 / Tf) * ((Kf / Tf) * Vf + Vr2) #16.12b
    output_ode[local_ix[4]] = (1.0 / Tr) * (V_th - Vm) #16.11

    #Update inner_vars
    get_inner_vars(device)[Vf_var] = Vf

    return
end

"""
Model of AVR Type II in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_avr_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V0_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, PSY.AVRTypeII)

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm = internal_states[4]

    #Define external states for device
    V_th = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)
    Vs = get_inner_vars(device)[V_pss_var]

    #Get parameters
    avr = PSY.get_avr(device)
    K0 = PSY.get_K0(avr)
    T1 = PSY.get_T1(avr)
    T2 = PSY.get_T2(avr)
    T3 = PSY.get_T3(avr)
    T4 = PSY.get_T4(avr)
    Te = PSY.get_Te(avr)
    Tr = PSY.get_Tr(avr)
    Vr_max = PSY.get_Vr_max(avr)
    Vr_min = PSY.get_Vr_min(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)

    #Compute auxiliary parameters
    Se_Vf = Ae * exp(Be * abs(Vf)) #16.13
    V_ref = V0_ref + Vs
    Vr = K0 * Vr2 + (T4 / T3) * (Vr1 + K0 * (T2 / T1) * (V_ref - Vm)) #16.21

    #Set anti-windup for Vr1. #TODO in callbacks
    #if Vr > Vr_max
    #    Vr = Vr_max
    #elseif Vr < Vr_min
    #    Vr = Vr_min
    #end

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = -(1.0 / Te) * (Vf * (1.0 + Se_Vf) - Vr) #16.18
    output_ode[local_ix[2]] = (1.0 / T1) * (K0 * (1.0 - (T2 / T1)) * (V_ref - Vm) - Vr1) #16.14
    output_ode[local_ix[3]] =
        (1.0 / (K0 * T3)) *
        ((1.0 - (T4 / T3)) * (Vr1 + K0 * (T2 / T1) * (V_ref - Vm)) - K0 * Vr2)  #16.20
    output_ode[local_ix[4]] = (1.0 / Tr) * (V_th - Vm) #16.11

    #Update inner_vars
    get_inner_vars(device)[Vf_var] = Vf

    return
end

"""
Model of Excitation System AC1A in Julia.
Refer to PSSE diagram block for details.
"""
function mdl_avr_ode!(
    device_states,
    output_ode,
    device::PSY.DynamicGenerator{M, S, PSY.AC1A, TG, P},
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}

    #Obtain references
    V_ref = PSY.get_ext(device)[CONTROL_REFS][V_ref_index]

    #Obtain avr
    avr = PSY.get_avr(device)

    #Obtain indices for component w/r to device
    local_ix = get_local_state_ix(device, typeof(avr))

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vm = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Ve = internal_states[4]
    Vr3 = internal_states[5]

    #Define external states for device
    V_th = sqrt(get_inner_vars(device)[VR_gen_var]^2 + get_inner_vars(device)[VI_gen_var]^2)
    Vs = get_inner_vars(device)[V_pss_var]
    Xad_Ifd = get_inner_vars(device)[Xad_Ifd_var]

    #Get parameters
    Tr = PSY.get_Tr(avr)
    Tb = PSY.get_Tb(avr)
    Tc = PSY.get_Tc(avr)
    Ka = PSY.get_Ka(avr)
    Ta = PSY.get_Ta(avr)
    #Va_min, Va_max = PSY.get_Va_lim(avr) #Not used without UEL or OEL
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

    #Compute auxiliary parameters
    I_N = Kc * Xad_Ifd / Ve
    V_FE = Kd * Xad_Ifd + Ke * Ve
    V_F = Vr3 + (Kf / Tf) * V_FE
    V_in = V_ref - Vm - V_F
    V_out = Vr1 + (Tc / Tb) * V_in
    Vf = Ve * rectifier_function(I_N)
    V_R = Vr2

    #Set anti-windup for Vr2. 
    if Vr2 > Vr_max
        V_R = Vr_max
    elseif Vr2 < Vr_min
        V_R = Vr_min
    end

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = (1.0 / Tr) * (V_th - Vm) #dVm/dt
    output_ode[local_ix[2]] = (1.0 / Tb) * (V_in * (1 - Tc / Tb) - Vr1) #dVr1/dt
    output_ode[local_ix[3]] = (1.0 / Ta) * (Ka * V_out - Vr2) #dVr2/dt
    output_ode[local_ix[4]] = (1.0 / Te) * (V_R - V_FE) #dVe/dt
    output_ode[local_ix[5]] = (1.0 / Tf) * (-(Kf / Tf) * V_FE - Vr3) #dVr3/dt

    #Update inner_vars
    get_inner_vars(device)[Vf_var] = Vf

    return
end
