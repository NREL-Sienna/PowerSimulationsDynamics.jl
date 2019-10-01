function mdl_avr_ode!(device_states,
                    output_ode,
                    device::DynGenerator{M, S, AVRFixed, TG, P})  where {M <: Machine,
                                                                         S <: Shaft,
                                                                         TG <: TurbineGov,
                                                                         P <: PSS}

    #Obtain parameters
    Vf = device.avr.Emf

    #Update inner vars
    device.inner_vars[Vf_var] = Vf

    return
end


function mdl_avr_ode!(device_states,
                      output_ode,
                      device::DynGenerator{M, S, AVRSimple, TG, P})  where {M <: Machine,
                                                                            S <: Shaft,
                                                                            TG <: TurbineGov,
                                                                            P <: PSS}

    #Obtain references
    V_ref = get_V_ref(device)

    local_ix = device.local_state_ix[device.avr]

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]

    #Define external states for device
    V_th = sqrt(device.inner_vars[VR_gen_var]^2 + device.inner_vars[VI_gen_var]^2)

    #Get Parameters
    Kv = device.avr.Kv

    #Compute ODEs
    output_ode[local_ix[1]] = Kv*(V_ref - V_th)

    #Update inner_vars
    device.inner_vars[Vf_var] = Vf

    return
end


"""
Model of AVR Type I in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_avr_ode!(device_states,
                      output_ode,
                      device::DynGenerator{M, S, AVRTypeI, TG, P})  where {M <: Machine,
                                                                           S <: Shaft,
                                                                           TG <: TurbineGov,
                                                                           P <: PSS}

    #Obtain references
    V0_ref = get_V_ref(device)

    local_ix = device.local_state_ix[device.avr]

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm  = internal_states[4]

    #Define external states for device
    V_th = sqrt(device.inner_vars[VR_gen_var]^2 + device.inner_vars[VI_gen_var]^2)
    Vs = device.inner_vars[V_pss_var]


    #Get parameters
    Ka = device.avr.Ka
    Ke = device.avr.Ke
    Kf = device.avr.Kf
    Ta = device.avr.Ta
    Te = device.avr.Te
    Tf = device.avr.Tf
    Tr = device.avr.Tr
    Vr_max = device.avr.Vr_max
    Vr_min = device.avr.Vr_min
    Ae = device.avr.Ae
    Be = device.avr.Be

    #Compute auxiliary parameters
    Se_Vf = Ae*exp(Be*abs(Vf)) #16.13
    V_ref = V0_ref + Vs

    #Set anti-windup for Vr1. NOT WORKING FOR INITIALIZATION USING NLSOLVE
    #if Vr1 > Vr_max
    #    Vr1 = Vr_max
    #elseif Vr1 < Vr_min
    #    Vr1 = Vr_min
    #end

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = -(1.0/Te)*(Vf*(Ke + Se_Vf) - Vr1) #16.12c
    output_ode[local_ix[2]] = (1.0/Ta)*(Ka*(V_ref - Vm - Vr2 - (Kf/Tf)*Vf) - Vr1) #16.12a
    output_ode[local_ix[3]] = -(1.0/Tf)*((Kf/Tf)*Vf + Vr2) #16.12b
    output_ode[local_ix[4]] = (1.0/Tr)*(V_th - Vm) #16.11

    #Update inner_vars
    device.inner_vars[Vf_var] = Vf

    return
end



"""
Model of AVR Type II in Julia.
Refer to Power System Modelling and Scripting by F. Milano for the equations
"""

function mdl_avr_ode!(device_states,
                      output_ode,
                      device::DynGenerator{M, S, AVRTypeII, TG, P})  where {M <: Machine,
                                                                           S <: Shaft,
                                                                           TG <: TurbineGov,
                                                                           P <: PSS}

    #Obtain references
    V0_ref = get_V_ref(device)

    local_ix = device.local_state_ix[device.avr]

    #Define inner states for component
    internal_states = @view device_states[local_ix]
    Vf = internal_states[1]
    Vr1 = internal_states[2]
    Vr2 = internal_states[3]
    Vm  = internal_states[4]

    #Define external states for device
    V_th = sqrt(device.inner_vars[VR_gen_var]^2 + device.inner_vars[VI_gen_var]^2)
    Vs = device.inner_vars[V_pss_var]


    #Get parameters
    K0 = device.avr.K0
    T1 = device.avr.T1
    T2 = device.avr.T2
    T3 = device.avr.T3
    T4 = device.avr.T4
    Te = device.avr.Te
    Tr = device.avr.Tr
    Vr_max = device.avr.Vr_max
    Vr_min = device.avr.Vr_min
    Ae = device.avr.Ae
    Be = device.avr.Be

    #Compute auxiliary parameters
    Se_Vf = Ae*exp(Be*abs(Vf)) #16.13
    V_ref = V0_ref + Vs
    Vr = K0*Vr2 + (T4/T3)*(Vr1 + K0*(T2/T1)*(V_ref - Vm)) #16.21

    #Set anti-windup for Vr1. NOT WORKING FOR INITIALIZATION USING NLSOLVE
    #if Vr > Vr_max
    #    Vr = Vr_max
    #elseif Vr < Vr_min
    #    Vr = Vr_min
    #end

    #Compute 4 States AVR ODE:
    output_ode[local_ix[1]] = -(1.0/Te)*(Vf*(1.0 + Se_Vf) - Vr) #16.18
    output_ode[local_ix[2]] = (1.0/T1)*(K0*(1.0 - (T2/T1))*(V_ref - Vm) - Vr1) #16.14
    output_ode[local_ix[3]] = (1.0/(K0*T3))*((1.0 - (T4/T3))*(Vr1 + K0*(T2/T1)*(V_ref - Vm)) - K0*Vr2)  #16.20
    output_ode[local_ix[4]] = (1.0/Tr)*(V_th - Vm) #16.11

    #Update inner_vars
    device.inner_vars[Vf_var] = Vf

    return
end
