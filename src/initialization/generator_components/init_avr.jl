function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = inner_vars[Vf_var]
    #Update Control Refs
    avr = PSY.get_avr(dynamic_device)
    set_V_ref(dynamic_device, Vf)
    PSY.set_Vf!(avr, Vf)
    PSY.set_V_ref!(avr, Vf)
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRSimple, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    #V_ref = Vm
    #Set Vf state equals to Vf0
    avr_ix = get_local_state_ix(dynamic_device, PSY.AVRSimple)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vf0 #δ
    #Set V_ref
    PSY.set_V_ref!(PSY.get_avr(dynamic_device), Vm)
    #Update Control Refs
    set_V_ref(dynamic_device, Vm)
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeI, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

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
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeI)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #δ
        avr_states[2] = sol_x0[2] #ω
        avr_states[3] = sol_x0[3]
        avr_states[4] = Vm
    end
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRTypeII, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    K0 = PSY.get_K0(avr)
    T1 = PSY.get_T1(avr)
    T2 = PSY.get_T2(avr)
    T3 = PSY.get_T3(avr)
    T4 = PSY.get_T4(avr)
    Te = PSY.get_Te(avr)
    Tr = PSY.get_Tr(avr)
    Ae = PSY.get_Ae(avr)
    Be = PSY.get_Be(avr)
    #Obtain saturated Vf
    Se_Vf = Ae * (exp(Be * abs(Vf0)))
    Va_min, Va_max = PSY.get_Va_lim(avr)

    #States of AVRTypeII are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr1, Vr2
    function f!(out, x)
        V_ref = x[1]
        Vr1 = x[2]
        Vr2 = x[3]

        y_ll1, dVr1_dt = lead_lag(V_ref - Vm, Vr1, K0, T2, T1)
        y_ll2, dVr2_dt = lead_lag(y_ll1, K0 * Vr2, 1.0, K0 * T4, K0 * T3)
        Vr = y_ll2
        _, dVf_dt = low_pass_modified(Vr, Vf0, 1.0, 1.0 + Se_Vf, Te)

        out[1] = dVf_dt #16.18
        out[2] = dVr1_dt #16.14
        out[3] = dVr2_dt #16.20
    end
    x0 = [1.0, Vf0, Vf0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.AVRTypeII)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #Vf
        avr_states[2] = sol_x0[2] #Vr1
        avr_states[3] = sol_x0[3] #Vr2
        avr_states[4] = Vm #Vm
        y_ll1, _ = lead_lag(sol_x0[1] - Vm, sol_x0[2], K0, T2, T1)
        y_ll2, _ = lead_lag(y_ll1, K0 * sol_x0[3], 1.0, K0 * T4, K0 * T3)
        if (y_ll2 > Va_max) || (y_ll2 < Va_min)
            @error("Regulator Voltage V_r = $(y_ll2) outside the limits")
        end
    end
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ESAC1A, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain I_fd obtained from Machine:
    Xad_Ifd0 = inner_vars[Xad_Ifd_var]
    #Obtain measured terminal voltage
    Vm0 = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

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
    sol = NLsolve.nlsolve(f_Ve!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
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
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
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

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.SEXS, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vm = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ta_Tb = PSY.get_Ta_Tb(avr)
    K = PSY.get_K(avr)
    V_min, V_max = PSY.get_V_lim(avr)

    #States of AVRTypeI are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr
    function f!(out, x)
        V_ref = x[1]
        Vr = x[2]

        V_in = V_ref - Vm
        V_LL = Vr + Ta_Tb * V_in

        out[1] = K * V_LL - Vf0
        out[2] = V_in * (1 - Ta_Tb) - Vr
    end
    x0 = [1.0, Vf0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        if (sol_x0[2] >= V_max + BOUNDS_TOLERANCE) ||
           (sol_x0[2] <= V_min - BOUNDS_TOLERANCE)
            @error(
                "Vr limits for AVR in $(PSY.get_name(dynamic_device)) (Vr = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min.  Consider updating the operating point."
            )
        end
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.SEXS)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #Vf
        avr_states[2] = sol_x0[2] #Vr
    end
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.SCRX, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain Ifd limiter
    Ifd = inner_vars[Xad_Ifd_var] # read Lad Ifd (field current times Lad)
    #Obtain measured terminal voltage
    Vm = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ta_Tb = PSY.get_Ta_Tb(avr)
    Tb = PSY.get_Tb(avr)
    Ta = Tb * Ta_Tb
    # dont need Te >> no derivative so block is 0 (multiply by K in steady state)
    K = PSY.get_K(avr)
    V_min, V_max = PSY.get_Efd_lim(avr) #Efd_lim (V_lim) **n
    switch = PSY.get_switch(avr) # reads switch parameters **n
    rc_rfd = PSY.get_rc_rfd(avr)

    # do the negative current and switch? or is that alr counted in? 
    #States of AVRTypeI are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr
    function f!(out, x)
        V_ref = x[1]
        Vr1 = x[2]

        V_in = V_ref - Vm # assume Vs is 0 when init
        #lead lag block
        #V_LL = Vr2 + Ta_Tb * V_in
        V_LL, dVr1_dt = lead_lag_mass_matrix(V_in, Vr1, 1.0, Ta, Tb) # 1st block
        Vr2 = K * V_LL

        # Switch multiplier
        mult = switch == 0 ? Vm : 1.0
        V_ex = mult * Vr2

        # negative current logic
        if rc_rfd == 0.0 # a float
            Efd = V_ex
        else
            if Ifd > 0.0
                Efd = V_ex
            else
                Efd = -Ifd * rc_rfd
            end
        end

        #V_ll output first block 
        out[1] = Efd - Vf0 # we are asking for Vf0 
        out[2] = dVr1_dt # make derivative 0 << tries to make it 0
        #out[2] = V_in * (1 - Ta_Tb) - Vr
    end # solve for Vref

    x0 = [1.0, Vf0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else # if converge
        sol_x0 = sol.zero
        Vr2_0 = (sol_x0[2] + Ta_Tb * (sol_x0[1] - Vm)) * K # K * V_LL
        #check the limits
        if (Vr2_0 >= V_max + BOUNDS_TOLERANCE) || (Vr2_0 <= V_min - BOUNDS_TOLERANCE)
            @error(
                "Vr limits for AVR in $(PSY.get_name(dynamic_device)) (Vr = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min.  Consider updating the operating point."
            )
        end
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.SCRX)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = sol_x0[2] #Vr1 // not Vf
        avr_states[2] = Vr2_0 #Vr2 // not Vr
    end
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.EXST1, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vt = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    #Obtain field winding current 
    Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base (for the available generator models)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Ka = PSY.get_Ka(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    Tc = PSY.get_Tc(avr)
    Tb = PSY.get_Tb(avr)
    Kc = PSY.get_Kc(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)

    # Check limits to field voltage 
    if (Vt * Vr_min - Kc * Ifd > Vf0) || (Vf0 > Vt * Vr_max - Kc * Ifd)
        @error(
            "Field Voltage for AVR in $(PSY.get_name(dynamic_device)) is $(Vf0) pu, which is outside its limits.  Consider updating the operating point."
        )
    end

    #Update V_ref
    Vref0 = Vt + Vf0 / Ka

    PSY.set_V_ref!(avr, Vref0)
    set_V_ref(dynamic_device, Vref0)

    #States of EXST1_PTI are Vm, Vll, Vr, Vfb

    #Update AVR states
    avr_ix = get_local_state_ix(dynamic_device, PSY.EXST1)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vt
    avr_states[2] = (1.0 - Tc / Tb) * Vf0 / Ka
    avr_states[3] = Vf0
    avr_states[4] = -Kf / Tf * Vf0
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.EXAC1, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vt0 = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    #Obtain field winding current
    Ifd0 = inner_vars[Xad_Ifd_var]

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Tb = PSY.get_Tb(avr)
    Tc = PSY.get_Tc(avr)
    Ka = PSY.get_Ka(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    Kc = PSY.get_Kc(avr)
    Kd = PSY.get_Kd(avr)
    Ke = PSY.get_Ke(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)

    #Solve Ve from rectifier function
    function f_Ve!(out, x)
        Ve = x[1]
        IN = Kc * Ifd0 / Ve

        out[1] = Vf0 - Ve * rectifier_function(IN)
    end
    x0 = [10.0] # initial guess for Ve
    sol = NLsolve.nlsolve(f_Ve!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else
        sol_x0 = sol.zero
        Ve = sol_x0[1]
        IN = Kc * Ifd0 / Ve
    end
    Se = saturation_function(avr, Ve) #To do saturation function
    VFE = Kd * Ifd0 + Ke * Ve + Se * Ve
    Vr2 = VFE
    if (Vr2 > Vr_max) || (Vr2 < Vr_min)
        @error("Regulator Voltage V_R = $(Vr2) outside the limits")
    end
    Vr3 = -(Kf / Tf) * VFE
    Tc_Tb_ratio = Tb <= eps() ? 0.0 : Tc / Tb
    Vr1 = (1 - Tc_Tb_ratio) * (VFE / Ka)
    Vm = Vt0
    Vref0 = Vt0 + VFE / Ka

    #Update V_ref
    PSY.set_V_ref!(avr, Vref0)
    set_V_ref(dynamic_device, Vref0)

    #States of EXAC1 are Vm, Vr1, Vr2, Ve, Vr3

    #Update AVR states
    avr_ix = get_local_state_ix(dynamic_device, typeof(avr))
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vm  #Vm
    avr_states[2] = Vr1 #Vr1
    avr_states[3] = Vr2 #Vr2
    avr_states[4] = Ve  #Ve
    avr_states[5] = Vr3 #Vr3
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ESST1A, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vt = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    #Obtain field winding current 
    Ifd = inner_vars[Xad_Ifd_var] # machine's field current in exciter base (for the available generator models)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Tc = PSY.get_Tc(avr)
    Tb = PSY.get_Tb(avr)
    Tc1 = PSY.get_Tc1(avr)
    Tb1 = PSY.get_Tb1(avr)
    Ka = PSY.get_Ka(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kc = PSY.get_Kc(avr)
    Kf = PSY.get_Kf(avr)
    Tf = PSY.get_Tf(avr)
    K_lr = PSY.get_K_lr(avr)
    I_lr = PSY.get_I_lr(avr)

    # Check limits to field voltage 
    if (Vt * Vr_min > Vf0) || (Vf0 > Vt * Vr_max - Kc * Ifd)
        @error(
            "Field Voltage for AVR in $(PSY.get_name(dynamic_device)) is $(Vf0) pu, which is outside its limits.  Consider updating the operating point."
        )
    end

    #Compute auxiliary parameters
    Itemp = K_lr * (Ifd - I_lr)
    Iresult = Itemp > 0 ? Itemp : 0

    Va = Vf0 + Iresult

    #Update V_ref
    Vref0 = Vt + Va / Ka

    PSY.set_V_ref!(avr, Vref0)
    set_V_ref(dynamic_device, Vref0)

    #States of ESST1A_PTI are Vm, Vr1, Vr2, Va, Vr3

    #Update AVR states
    avr_ix = get_local_state_ix(dynamic_device, PSY.ESST1A)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vt
    avr_states[2] = (1 - Tc / Tb) * Va / Ka
    avr_states[3] = (1 - Tc1 / Tb1) * Va / Ka
    avr_states[4] = Va
    avr_states[5] = -Kf / Tf * Vf0
    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ST6B, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain measured terminal voltage
    Vt0 = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)
    #Obtain field winding current
    Ifd0 = inner_vars[Xad_Ifd_var]

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    Tr = PSY.get_Tr(avr)
    K_pa = PSY.get_K_pa(avr) #k_pa>0
    K_ia = PSY.get_K_ia(avr)
    K_da = PSY.get_K_da(avr)
    T_da = PSY.get_T_da(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    K_ff = PSY.get_K_ff(avr)
    K_m = PSY.get_K_m(avr)
    K_ci = PSY.get_K_ci(avr) #K_cl in pss
    K_lr = PSY.get_K_lr(avr)
    I_lr = PSY.get_I_lr(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kg = PSY.get_Kg(avr)
    Tg = PSY.get_Tg(avr) #T_g>0

    #Compute auxiliary parameters
    Vr = Vf0 / Vt0
    Vg = Kg * Vf0

    Vr2 = max(Vr_min, ((I_lr * K_ci) - Ifd0) * K_lr)
    if (Vr - Vr2 > eps()) # or using STRICT_NLSOLVE_F_TOLERANCE?

        #Check saturation

        Va = (Vr + K_m * Vg) / (K_ff + K_m)

        #Check ss error according to control parameters
        if K_ia < eps()
            Vref0 = Va / K_pa + Vt0
            error(
                "AVR in $(PSY.get_name(dynamic_device)) has non positive integrator gain. Please update it.",
            )
        else
            Vref0 = Vt0
        end

        if !((Vr < Vr_max) && (Vr > Vr_min))
            @error(
                "AVR controller in $(PSY.get_name(dynamic_device)) is saturated. Consider updating the operating point."
            )
        end
    else
        error(
            "Current limiter of AVR in $(PSY.get_name(dynamic_device)) is activated. Consider updating the operating point.",
        )
    end

    PSY.set_V_ref!(avr, Vref0)
    set_V_ref(dynamic_device, Vref0)
    @warn(
        "I_LR parameter was updated from $(I_lr) to $(Ifd0) of $(PSY.get_name(dynamic_device)) to have zero current field error."
    )
    PSY.set_I_lr!(avr, Ifd0 * 10.0)

    #States of ST6B [:Vm, :x_i, :x_d, :Vg]
    #Update AVR states
    avr_ix = get_local_state_ix(dynamic_device, PSY.ST6B)
    avr_states = @view device_states[avr_ix]
    avr_states[1] = Vt0
    avr_states[2] = Va / K_ia
    avr_states[3] = 0
    avr_states[4] = Vg

    return
end

function initialize_avr!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.ST8C, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #Obtain Vf0 solved from Machine
    Vf0 = inner_vars[Vf_var]
    #Obtain Ifd limiter
    Ifd = inner_vars[Xad_Ifd_var] # read Lad Ifd (field current times Lad)
    #Obtain measured terminal voltage
    Vt0 = sqrt(inner_vars[VR_gen_var]^2 + inner_vars[VI_gen_var]^2)

    #Get parameters
    avr = PSY.get_avr(dynamic_device)
    SW1_Flag = PSY.get_SW1_Flag(avr)
    if SW1_Flag == 1
        error("Source from generator terminal voltage not supported.")
    end
    Tr = PSY.get_Tr(avr)
    K_pr = PSY.get_K_pr(avr)
    K_ir = PSY.get_K_ir(avr)
    Vpi_min, Vpi_max = PSY.get_Vpi_lim(avr)
    K_pa = PSY.get_K_pa(avr)
    K_ia = PSY.get_K_ia(avr)
    Va_min, Va_max = PSY.get_Va_lim(avr)
    K_a = PSY.get_K_a(avr)
    T_a = PSY.get_T_a(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    K_f = PSY.get_K_f(avr)
    T_f = PSY.get_T_f(avr)
    K_c1 = PSY.get_K_c1(avr)
    K_p = PSY.get_K_p(avr)
    K_i2 = PSY.get_K_i2(avr)
    VB1_max = PSY.get_VB1_max(avr)

    if K_i2 != 0.0
        error("Feedforward Current for AVR ST8C not implemented yet.")
    end
    if K_ia == 0.0
        error(
            "Integrator gain cannot be zero for AVR ST8C of $(PSY.get_name(dynamic_device))",
        )
    end
    #To solve V_ref, Vr
    function f!(out, x)
        V_ref = x[1]
        Vm = x[2] # Sensed Voltage
        x_a1 = x[3] # Regulator Integrator
        x_a2 = x[4] # Field Regulator
        x_a3 = x[5] # Controller Integrator
        x_a4 = x[6] # Regulator Feedback

        #TODO: Implement Terminal Current FF for AVR
        V_b2 = 0.0
        #TODO: Implement Voltage Compensation if needed
        V_e = K_p

        # Compute block derivatives
        Vs = 0.0
        V_pi_in = V_ref + Vs - Vm
        Ifd_ref, dxa1_dt = pi_block_nonwindup(V_pi_in, x_a1, K_pr, K_ir, Vpi_min, Vpi_max)
        Ifd_diff = Ifd_ref - x_a4
        pi_out, dxa2_dt = pi_block_nonwindup(Ifd_diff, x_a2, K_pa, K_ia, Va_min, Va_max)
        _, dxa3_dt = low_pass_nonwindup_mass_matrix(pi_out, x_a3, K_a, T_a, Vr_min, Vr_max)
        _, dxa4_dt = low_pass_mass_matrix(Ifd, x_a4, K_f, T_f)

        # Compute V_b1
        I_n1 = K_c1 * Ifd / V_e
        F_ex = rectifier_function(I_n1)
        V_b1 = F_ex * V_e

        Efd = V_b1 * x_a3 + V_b2

        #V_ll output first block 
        out[1] = Efd - Vf0 # we are asking for Vf0 
        out[2] = Vm - Vt0 # Low Pass Sensor
        out[3] = dxa1_dt
        out[4] = dxa2_dt
        out[5] = dxa3_dt
        out[6] = dxa4_dt
    end
    # Initial Guess
    V_e0 = K_p
    I_n10 = K_c1 * Ifd / V_e0
    F_ex0 = rectifier_function(I_n10)
    V_b10 = F_ex0 * V_e0
    x_a30 = Vf0 / V_b10
    x_a20 = x_a30 / (K_a * K_ia)
    x0 = [Vt0, Vt0, K_f * Ifd / K_ir, x_a20, x_a30, Ifd * K_f]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization of AVR in $(PSY.get_name(static)) failed")
    else # if converge
        sol_x0 = sol.zero
        Vreg = sol_x0[6]
        if (Vreg > Vr_max) || (Vreg < Vr_min)
            @error(
                "Regulator Voltage V_R = $(Vreg) outside the limits for AVR of $(PSY.get_name(dynamic_device)). Consider updating its limits."
            )
        end
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref(dynamic_device, sol_x0[1])
        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.ST8C)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = sol_x0[2]
        avr_states[2] = sol_x0[3]
        avr_states[3] = sol_x0[4]
        avr_states[4] = sol_x0[5]
        avr_states[5] = sol_x0[6]
    end
end
