function initialize_shaft!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.SingleMass, A, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS} end

function initialize_shaft!(
    device_states,
    device_parameters,
    static::PSY.StaticInjection,
    dynamic_device::DynamicWrapper{PSY.DynamicGenerator{M, PSY.FiveMassShaft, A, TG, P}},
    inner_vars::AbstractVector,
) where {M <: PSY.Machine, A <: PSY.AVR, TG <: PSY.TurbineGov, P <: PSY.PSS}
    shaft_ix = get_local_state_ix(dynamic_device, PSY.FiveMassShaft)
    shaft_states = @view device_states[shaft_ix]
    δ0 = shaft_states[1]
    δ = δ0
    ω = shaft_states[2]
    ω_hp = ω
    ω_ip = ω
    ω_lp = ω
    ω_ex = ω
    ω_sys = ω

    #Get parameters
    local_ix_params = get_local_parameter_ix(dynamic_device, PSY.FiveMassShaft)
    internal_params = @view device_parameters[local_ix_params]
    _,
    _,
    _,
    _,
    _,
    D,
    D_hp,
    D_ip,
    D_lp,
    D_ex,
    D_12,
    D_23,
    D_34,
    D_45,
    K_hp,
    K_ip,
    K_lp,
    K_ex = internal_params

    #Obtain inner vars
    τe = inner_vars[τe_var]
    τm0 = inner_vars[τm_var]

    function f!(out, x)
        τm = x[1]
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

    x0 = [τm0, δ0, δ0, δ0, δ0]
    sol = NLsolve.nlsolve(f!, x0; ftol = STRICT_NLSOLVE_F_TOLERANCE)
    if !NLsolve.converged(sol)
        @warn("Initialization in Shaft failed")
    else
        sol_x0 = sol.zero
        inner_vars[τm_var] = sol_x0[1] #τm
        shaft_states[3] = sol_x0[2] #δ_hp
        shaft_states[4] = ω #ω_hp
        shaft_states[5] = sol_x0[3] #δ_ip
        shaft_states[6] = ω #ω_ip
        shaft_states[7] = sol_x0[4] #δ_lp
        shaft_states[8] = ω #ω_lp
        shaft_states[9] = sol_x0[5] #δ_ex
        shaft_states[10] = ω #ω_ex
    end
    return
end
