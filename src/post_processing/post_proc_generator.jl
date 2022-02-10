"""
Function to obtain the output current time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64}
) where {G <: PSY.DynamicGenerator}

    #Obtain Data
    sys = get_system(res)

    #Get machine
    machine = PSY.get_machine(dynamic_device)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase
    return _machine_current(
        machine,
        PSY.get_name(dynamic_device),
        V_R,
        V_I,
        base_power_ratio,
        res,
        dt
    )
end

"""
Function to obtain the output current time series of a Classic Machine model out of the DAE Solution. It is dispatched via the machine type.

"""
function _machine_current(
    machine::PSY.BaseMachine,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)

    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    eq_p = PSY.get_eq_p(machine)

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ)
    I_I = similar(δ)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]
        #Obtain electric current
        i_dq[1] = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * (eq_p .- V_q) - R * V_d)  #15.36
        i_dq[2] = (1.0 / (R^2 + Xd_p^2)) * (Xd_p * V_d + R * (eq_p .- V_q)) #15.36

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a One-D-One-Q model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::PSY.OneDOneQMachine,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ed_p = post_proc_state_series(res, (name, :ed_p), dt)

    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xq_p = PSY.get_Xq_p(machine)

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]
        #Obtain electric current
        i_dq[1] =
            (1.0 / (R^2 + Xd_p * Xq_p)) * (Xq_p * (eq_p[ix] - V_q) + R * (ed_p[ix] - V_d))  #15.32
        i_dq[2] =
            (1.0 / (R^2 + Xd_p * Xq_p)) * (-Xd_p * (ed_p[ix] - V_d) + R * (eq_p[ix] - V_q))  #15.32

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a SimpleMarconato or SimpleAndersonFouad model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::PSY.Union{PSY.SimpleMarconatoMachine, PSY.SimpleAFMachine},
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_pp = post_proc_state_series(res, (name, :eq_pp), dt)
    _, ed_pp = post_proc_state_series(res, (name, :ed_pp), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]
        #Obtain electric current
        i_dq[1] =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (Xq_pp * (eq_pp[ix] - V_q) + R * (ed_pp[ix] - V_d))      #15.25
        i_dq[2] =
            (1.0 / (R^2 + Xd_pp * Xq_pp)) *
            (-Xd_pp * (ed_pp[ix] - V_d) + R * (eq_pp[ix] - V_q))      #15.25

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a Marconato or AndersonFouad model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::Union{PSY.MarconatoMachine, PSY.AndersonFouadMachine},
    name::String,
    ::Vector{Float64},
    ::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_pp = post_proc_state_series(res, (name, :eq_pp), dt)
    _, ed_pp = post_proc_state_series(res, (name, :ed_pp), dt)
    _, ψd = post_proc_state_series(res, (name, :ψd), dt)
    _, ψq = post_proc_state_series(res, (name, :ψq), dt)

    #Get parameters
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]

        #Obtain electric current
        i_dq[1] = (1.0 / Xd_pp) * (eq_pp[ix] - ψd[ix])      #15.18
        i_dq[2] = (1.0 / Xq_pp) * (-ed_pp[ix] - ψq[ix])     #15.18

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a GENROU/GENROE model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::Union{PSY.RoundRotorQuadratic, PSY.RoundRotorExponential},
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ed_p = post_proc_state_series(res, (name, :ed_p), dt)
    _, ψ_kd = post_proc_state_series(res, (name, :ψ_kd), dt)
    _, ψ_kq = post_proc_state_series(res, (name, :ψ_kq), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    ψq_pp = γ_q1 * ed_p + ψ_kq * (1 - γ_q1)
    ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψ_kd

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]

        #Obtain electric current
        i_dq[1] =
            (1.0 / (R^2 + Xq_pp * Xd_pp)) *
            (-R * (V_d - ψq_pp[ix]) + Xq_pp * (-V_q + ψd_pp[ix]))
        i_dq[2] =
            (1.0 / (R^2 + Xq_pp * Xd_pp)) *
            (Xd_pp * (V_d - ψq_pp[ix]) + R * (-V_q + ψd_pp[ix]))

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end

"""
Function to obtain the output current time series of a GENSAL/GENSAE model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::Union{PSY.SalientPoleQuadratic, PSY.SalientPoleExponential},
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64}
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ψ_kd = post_proc_state_series(res, (name, :ψ_kd), dt)
    _, ψq_pp = post_proc_state_series(res, (name, :ψq_pp), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)

    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]

        #Obtain electric current
        i_dq[1] =
            (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d + ψq_pp[ix]) + Xd_pp * (ψd_pp[ix] - V_q))
        i_dq[2] =
            (1.0 / (R^2 + Xd_pp^2)) * (Xd_pp * (V_d + ψq_pp[ix]) + R * (ψd_pp[ix] - V_q))

        I_R[ix], I_I[ix] = base_power_ratio * dq_ri(v) * i_dq
    end
    return ts, I_R, I_I
end
