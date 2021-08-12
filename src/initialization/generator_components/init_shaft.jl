function initialize_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.SingleMass, A, TG, P}},
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS} end

function initialize_shaft!(
    device_states,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.FiveMassShaft, A, TG, P}},
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    shaft_ix = get_local_state_ix(dynamic_device, PSY.FiveMassShaft)
    shaft_states = @view device_states[shaft_ix]
    δ0 = shaft_states[1]
    ω = shaft_states[2]
    ω_hp = ω
    ω_ip = ω
    ω_lp = ω
    ω_ex = ω
    ω_sys = ω

    #Obtain parameters
    shaft = PSY.get_shaft(dynamic_device)
    D = PSY.get_D(shaft)
    D_hp = PSY.get_D_hp(shaft)
    D_ip = PSY.get_D_ip(shaft)
    D_lp = PSY.get_D_lp(shaft)
    D_ex = PSY.get_D_ex(shaft)
    D_12 = PSY.get_D_12(shaft)
    D_23 = PSY.get_D_23(shaft)
    D_34 = PSY.get_D_34(shaft)
    D_45 = PSY.get_D_45(shaft)
    K_hp = PSY.get_K_hp(shaft)
    K_ip = PSY.get_K_ip(shaft)
    K_lp = PSY.get_K_lp(shaft)
    K_ex = PSY.get_K_ex(shaft)

    #Obtain inner vars
    τe = get_inner_vars(dynamic_device)[τe_var]
    τm = get_inner_vars(dynamic_device)[τm_var]

    function f!(out, x)
        δ = x[1]
        δ_hp = x[2]
        δ_ip = x[3]
        δ_lp = x[4]
        δ_ex = x[5]

        out[1] = (
            -τe - D * (ω - ω_sys) - D_34 * (ω - ω_lp) - D_45 * (ω - ω_ex) +
            K_lp * (δ_lp - δ) +
            K_ex * (δ_ex - δ)
        )

        out[2] = τm - D_hp * (ω_hp - ω_sys) - D_12 * (ω_hp - ω_ip) + K_hp * (δ_ip - δ_hp)

        out[3] = (
            -D_ip * (ω_ip - ω_sys) - D_12 * (ω_ip - ω_hp) - D_23 * (ω_ip - ω_lp) +
            K_hp * (δ_hp - δ_ip) +
            K_ip * (δ_lp - δ_ip)
        )

        out[4] = (
            -D_lp * (ω_lp - ω_sys) - D_23 * (ω_lp - ω_ip) - D_34 * (ω_lp - ω) +
            K_ip * (δ_ip - δ_lp) +
            K_lp * (δ - δ_lp)
        )
        out[5] = -D_ex * (ω_ex - ω_sys) - D_45 * (ω_ex - ω) + K_ex * (δ - δ_ex)
    end

    x0 = [δ0, δ0, δ0, δ0, δ0]
    sol = NLsolve.nlsolve(f!, x0)
    if !NLsolve.converged(sol)
        @warn("Initialization in Shaft failed")
    else
        sol_x0 = sol.zero
        shaft_states[1] = sol_x0[1] #δ0
        shaft_states[2] = ω
        shaft_states[3] = sol_x0[2] #δ_hp
        shaft_states[4] = ω #ω_hp
        shaft_states[5] = sol_x0[3] #δ_ip
        shaft_states[6] = ω #ω_ip
        shaft_states[7] = sol_x0[4] #δ_lp
        shaft_states[8] = ω #ω_lp
        shaft_states[9] = sol_x0[5] #δ_ex
        shaft_states[10] = ω #ω_ex
    end
end
