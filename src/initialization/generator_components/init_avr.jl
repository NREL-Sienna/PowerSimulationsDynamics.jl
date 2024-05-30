function initialize_avr!(
    device_states,
    p,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, S, PSY.AVRFixed, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    #In AVRFixed, V_ref is used as Vf
    Vf = inner_vars[Vf_var]
    #Update Control Refs
    avr = PSY.get_avr(dynamic_device)
    set_V_ref!(p, Vf)
    PSY.set_Vf!(avr, Vf)
    PSY.set_V_ref!(avr, Vf)
    return
end

function initialize_avr!(
    device_states,
    device_parameters,
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
    device_parameters,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.AVRTypeI)
    internal_params = @view device_parameters[local_ix_params]
    Ka, Ke, Kf, _, _, Tf, _, Ae, Be = internal_params
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
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
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
    device_parameters,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.AVRTypeII)
    internal_params = @view device_parameters[local_ix_params]
    K0,
    T1,
    T2,
    T3,
    T4,
    Te,
    _,
    Va_min,
    Va_max,
    Ae,
    Be = internal_params
    #Obtain saturated Vf
    Se_Vf = Ae * (exp(Be * abs(Vf0)))

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
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
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
            CRC.@ignore_derivatives @error(
                "Regulator Voltage V_r = $(y_ll2) outside the limits"
            )
        end
    end
    return
end

function initialize_avr!(
    device_states,
    device_parameters,
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
    #Get parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.ESAC1A)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Tb,
    Tc,
    Ka,
    Ta,
    Va_min,
    Va_max,
    Te,
    Kf,
    Tf,
    Kc,
    Kd,
    Ke,
    Vr_min,
    Vr_max = internal_params
    inv_Tr = Tr < eps() ? 1.0 : 1.0 / Tr
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
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        V_e0 = sol_x0[1]
        I_N0 = Kc * Xad_Ifd0 / V_e0
    end
    Se0 = saturation_function(avr, V_e0) #To do saturation function
    V_FE0 = Kd * Xad_Ifd0 + Ke * V_e0 + Se0 * V_e0
    V_r20 = V_FE0
    if (V_r20 > Vr_max) || (V_r20 < Vr_min)
        CRC.@ignore_derivatives @error(
            "Regulator Voltage V_R = $(V_r20) outside the limits"
        )
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
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
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
    p,
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
    params = p[:params][:AVR]
    V_min = params[:V_lim][:min]
    V_max = params[:V_lim][:max]

    #States of AVRTypeI are Vf, Vr1, Vr2, Vm
    #To solve V_ref, Vr
    function f!(out, x, params)
        V_ref = x[1]
        Vr = x[2]
        Ta_Tb = params[:Ta_Tb]
        K = params[:K]

        V_in = V_ref - Vm
        V_LL = Vr + Ta_Tb * V_in

        out[1] = K * V_LL - Vf0
        out[2] = V_in * (1 - Ta_Tb) - Vr
    end
    x0 = [1.0, Vf0]
    prob = NonlinearSolve.NonlinearProblem(f!, x0, params)
    sol = NonlinearSolve.solve(
        prob,
        NonlinearSolve.TrustRegion();
        reltol = STRICT_NLSOLVE_F_TOLERANCE,
        abstol = STRICT_NLSOLVE_F_TOLERANCE,
    )
    if !SciMLBase.successful_retcode(sol)
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.u
        if (sol_x0[2] >= V_max + BOUNDS_TOLERANCE) ||
           (sol_x0[2] <= V_min - BOUNDS_TOLERANCE)
            CRC.@ignore_derivatives @error(
                "Vr limits for AVR in $(PSY.get_name(dynamic_device)) (Vr = $(sol_x0[2])), outside its limits V_max = $V_max, Vmin = $V_min.  Consider updating the operating point."
            )
        end
        #Update V_ref
        PSY.set_V_ref!(avr, sol_x0[1])
        set_V_ref!(p, sol_x0[1])

        #Update AVR states
        avr_ix = get_local_state_ix(dynamic_device, PSY.SEXS)
        avr_states = @view device_states[avr_ix]
        avr_states[1] = Vf0 #Vf
        avr_states[2] = sol_x0[2] #Vr
    end
end

function initialize_avr!(
    device_states,
    device_parameters,
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
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
    else # if converge
        sol_x0 = sol.zero
        Vr2_0 = (sol_x0[2] + Ta_Tb * (sol_x0[1] - Vm)) * K # K * V_LL
        #check the limits
        if (Vr2_0 >= V_max + BOUNDS_TOLERANCE) || (Vr2_0 <= V_min - BOUNDS_TOLERANCE)
            CRC.@ignore_derivatives @error(
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
    device_parameters,
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
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.EXST1)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Vi_min,
    Vi_max,
    Tc,
    Tb,
    Ka,
    Ta,
    Vr_min,
    Vr_max,
    Kc,
    Kf,
    Tf = internal_params

    # Check limits to field voltage 
    if (Vt * Vr_min - Kc * Ifd > Vf0) || (Vf0 > Vt * Vr_max - Kc * Ifd)
        CRC.@ignore_derivatives @error(
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
    device_parameters,
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
    #Get parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.EXAC1)
    internal_params = @view device_parameters[local_ix_params]
    Tr,
    Tb,
    Tc,
    Ka,
    Ta,
    Vr_min,
    Vr_max,
    Te,
    Kf,
    Tf,
    Kc,
    Kd,
    Ke = internal_params

    #Solve Ve from rectifier function
    function f_Ve!(out, x)
        Ve = x[1]
        IN = Kc * Ifd0 / Ve

        out[1] = Vf0 - Ve * rectifier_function(IN)
    end
    x0 = [10.0] # initial guess for Ve
    sol = NLsolve.nlsolve(f_Ve!, x0)
    if !NLsolve.converged(sol)
        CRC.@ignore_derivatives @warn(
            "Initialization of AVR in $(PSY.get_name(static)) failed"
        )
    else
        sol_x0 = sol.zero
        Ve = sol_x0[1]
        IN = Kc * Ifd0 / Ve
    end
    Se = saturation_function(avr, Ve) #To do saturation function
    VFE = Kd * Ifd0 + Ke * Ve + Se * Ve
    Vr2 = VFE
    if (Vr2 > Vr_max) || (Vr2 < Vr_min)
        CRC.@ignore_derivatives @error("Regulator Voltage V_R = $(Vr2) outside the limits")
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
    device_parameters,
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
        CRC.@ignore_derivatives @error(
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
