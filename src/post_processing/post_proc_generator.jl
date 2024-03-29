"""
Default function to compute output current. 
Returns an error

"""
function compute_output_current(
    ::SimulationResults,
    dynamic_device::I,
    ::Vector{Float64},
    ::Vector{Float64},
    ::Union{Nothing, Float64, Vector{Float64}},
) where {I <: PSY.DynamicInjection}

    #Return error
    error(
        "Output current for device type $(typeof(dynamic_device)) is not implemented yet.",
    )
end

"""
Function to obtain the output current time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64, Vector{Float64}},
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
        dt,
    )
end

"""
Function to obtain the output current time series of a AggregateDistributedGenerationA (DERA) model out of the DAE Solution.
It receives the simulation inputs, the dynamic device and bus voltage.

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::PSY.AggregateDistributedGenerationA,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64, Vector{Float64}},
)

    #Obtain Data
    sys = get_system(res)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase
    Freq_Flag = PSY.get_Freq_Flag(dynamic_device)
    name = PSY.get_name(dynamic_device)
    if Freq_Flag == 1
        _, Pord = post_proc_state_series(res, (name, :Pord), dt)
        _, dPord = post_proc_state_series(res, (name, :dPord), dt)
        Tpord = PSY.get_Tpord(dynamic_device)
        P_lim = PSY.get_P_lim(dynamic_device)
    end

    # Get states
    θ = atan.(V_I ./ V_R)
    ts, Ip = post_proc_state_series(res, (name, :Ip), dt)
    ts, Iq = post_proc_state_series(res, (name, :Iq), dt)
    _, Mult = post_proc_state_series(res, (name, :Mult), dt)
    _, Vmeas = post_proc_state_series(res, (name, :Vmeas), dt)
    Iq_neg = -Iq
    Ip_cmd = Ip
    Iq_cmd = Iq

    # Get Params
    Tg = PSY.get_Tg(dynamic_device)
    rrpwr = PSY.get_rrpwr(dynamic_device)
    P_ref = PSY.get_P_ref(dynamic_device)

    I_R = similar(Ip)
    I_I = similar(Iq)
    for (ix, Ip_cmd_val) in enumerate(Ip_cmd)
        Ip_min, Ip_max, _, _ = current_limit_logic(dynamic_device, Ip_cmd_val, Iq_cmd[ix])
        if Ip[ix] >= 0
            Rup = abs(rrpwr)
            Rdown = -Inf
        else
            Rdown = -abs(rrpwr)
            Rup = Inf
        end
        if Freq_Flag == 0
            Ip_input = clamp(P_ref / max(Vmeas[ix], 0.01), Ip_min, Ip_max) * Mult[ix]
            Ip_limited, _ = low_pass_nonwindup_ramp_limits(
                Ip_input,
                Ip[ix],
                1.0,
                Tg,
                -Inf,
                Inf,
                Rdown,
                Rup,
            )
        else
            Pord_limited, _ = low_pass_nonwindup_mass_matrix(
                dPord[ix],
                Pord[ix],
                1.0,
                Tpord,
                P_lim[:min],
                P_lim[:max],
            )
            Ip_input = clamp(Pord_limited / max(Vmeas[ix], 0.01), Ip_min, Ip_max) * Mult[ix]
            Ip_limited, _ = low_pass_nonwindup_ramp_limits(
                Ip_input,
                Ip[ix],
                1.0,
                Tg,
                -Inf,
                Inf,
                Rdown,
                Rup,
            )
        end
        I_R[ix] = real(complex(Ip_limited, Iq_neg[ix]) * exp(im * θ[ix])) * base_power_ratio
        I_I[ix] = imag(complex(Ip_limited, Iq_neg[ix]) * exp(im * θ[ix])) * base_power_ratio
    end

    return ts, I_R, I_I
end

"""
Function to obtain the field current time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_field_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {G <: PSY.DynamicGenerator}

    #Obtain Data
    sys = get_system(res)

    #Get machine
    machine = PSY.get_machine(dynamic_device)
    return _field_current(machine, PSY.get_name(dynamic_device), V_R, V_I, res, dt)
end

"""
Function to obtain the field voltage time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific voltage.

"""
function compute_field_voltage(
    res::SimulationResults,
    dynamic_device::G,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {G <: PSY.DynamicGenerator}

    #Get AVR
    avr = PSY.get_avr(dynamic_device)
    return _field_voltage(avr, PSY.get_name(dynamic_device), res, dt)
end

"""
Function to obtain the pss output time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific output.

"""
function compute_pss_output(
    res::SimulationResults,
    dynamic_device::G,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {G <: PSY.DynamicGenerator}

    #Get PSS
    pss = PSY.get_pss(dynamic_device)
    return _pss_output(pss, PSY.get_name(dynamic_device), res, dt)
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific torque.

"""
function compute_mechanical_torque(
    res::SimulationResults,
    dynamic_device::G,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {G <: PSY.DynamicGenerator}

    #Get TG
    tg = PSY.get_prime_mover(dynamic_device)
    return _mechanical_torque(tg, PSY.get_name(dynamic_device), res, dt)
end

"""
Function to obtain the output frequency time series of a DynamicGenerator

"""
function compute_frequency(
    res::SimulationResults,
    dyn_device::G,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {G <: PSY.DynamicGenerator}
    name = PSY.get_name(dyn_device)
    ts, ω = post_proc_state_series(res, (name, :ω), dt)
    return ts, ω
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
    dt::Union{Nothing, Float64, Vector{Float64}},
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
    dt::Union{Nothing, Float64, Vector{Float64}},
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
    dt::Union{Nothing, Float64, Vector{Float64}},
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
    dt::Union{Nothing, Float64, Vector{Float64}},
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
Function to obtain the output current time series of a SauerPaiMachine model out of the DAE Solution. It is dispatched via the machine type.
"""
function _machine_current(
    machine::PSY.SauerPaiMachine,
    name::String,
    ::Vector{Float64},
    ::Vector{Float64},
    base_power_ratio::Float64,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ed_p = post_proc_state_series(res, (name, :ed_p), dt)
    _, ψd = post_proc_state_series(res, (name, :ψd), dt)
    _, ψq = post_proc_state_series(res, (name, :ψq), dt)
    _, ψd_pp = post_proc_state_series(res, (name, :ψd_pp), dt)
    _, ψq_pp = post_proc_state_series(res, (name, :ψq_pp), dt)

    #Get parameters
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = PSY.get_Xq_pp(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)

    i_dq = Vector{Float64}(undef, 2)
    I_R = similar(δ, Float64)
    I_I = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]

        #Obtain electric current
        i_dq[1] = (1.0 / Xd_pp) * (γ_d1 * eq_p[ix] - ψd[ix] + (1 - γ_d1) * ψd_pp[ix])      #15.15
        i_dq[2] = (1.0 / Xq_pp) * (-γ_q1 * ed_p[ix] - ψq[ix] + (1 - γ_q1) * ψq_pp[ix])     #15.15

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
    dt::Union{Nothing, Float64, Vector{Float64}},
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
    dt::Union{Nothing, Float64, Vector{Float64}},
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

"""
Function to obtain the field current time series of a Dynamic Generator. It is dispatched via the machine type.
By default, machine does not have support for field current

"""
function _field_current(
    machine::M,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {M <: PSY.Machine}
    @warn("Field current is not supported in the machine type $(M). Returning zeros.")
    ts, _ = post_proc_state_series(res, (name, :δ), dt)
    I_fd = zeros(length(ts))
    return ts, I_fd
end

"""
Function to obtain the field current time series of a Dynamic Generator with machine type GENROU/GENROE.

"""
function _field_current(
    machine::M,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {M <: Union{PSY.RoundRotorQuadratic, PSY.RoundRotorExponential}}
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ed_p = post_proc_state_series(res, (name, :ed_p), dt)
    _, ψ_kd = post_proc_state_series(res, (name, :ψ_kd), dt)
    _, ψ_kq = post_proc_state_series(res, (name, :ψ_kd), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #Additional Fluxes
    ψq_pp = γ_q1 * ed_p + ψ_kq * (1 - γ_q1)
    ψd_pp = γ_d1 * eq_p + γ_d2 * (Xd_p - Xl) * ψ_kd
    ψ_pp = sqrt.(ψd_pp .^ 2 + ψq_pp .^ 2)

    I_fd = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]

        #Obtain electric current
        I_d =
            (1.0 / (R^2 + Xq_pp * Xd_pp)) *
            (-R * (V_d - ψq_pp[ix]) + Xq_pp * (-V_q + ψd_pp[ix]))

        Se = saturation_function(machine, ψ_pp[ix])
        I_fd[ix] =
            eq_p[ix] +
            (Xd - Xd_p) * (γ_d1 * I_d - γ_d2 * ψ_kd[ix] + γ_d2 * eq_p[ix]) +
            Se * ψd_pp[ix]
    end
    return ts, I_fd
end

"""
Function to obtain the field current time series of a Dynamic Generator with machine type GENSAL.

"""
function _field_current(
    machine::PSY.SalientPoleQuadratic,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ψ_kd = post_proc_state_series(res, (name, :ψ_kd), dt)
    _, ψq_pp = post_proc_state_series(res, (name, :ψq_pp), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd

    I_fd = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]

        #Obtain electric current
        I_d = (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d + ψq_pp[ix]) + Xd_pp * (ψd_pp[ix] - V_q))
        Se = saturation_function(machine, eq_p[ix])
        I_fd[ix] =
            eq_p[ix] +
            Se * eq_p[ix] +
            (Xd - Xd_p) * (I_d + γ_d2 * (eq_p[ix] - ψ_kd[ix] - (Xd_p - Xl) * I_d))
    end
    return ts, I_fd
end

"""
Function to obtain the field current time series of a Dynamic Generator with machine type GENSAE.

"""
function _field_current(
    machine::PSY.SalientPoleExponential,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, δ = post_proc_state_series(res, (name, :δ), dt)
    _, eq_p = post_proc_state_series(res, (name, :eq_p), dt)
    _, ψ_kd = post_proc_state_series(res, (name, :ψ_kd), dt)
    _, ψq_pp = post_proc_state_series(res, (name, :ψq_pp), dt)

    #Get parameters
    R = PSY.get_R(machine)
    Xd = PSY.get_Xd(machine)
    Xd_p = PSY.get_Xd_p(machine)
    Xd_pp = PSY.get_Xd_pp(machine)
    Xq_pp = Xd_pp
    Xl = PSY.get_Xl(machine)
    γ_d1 = PSY.get_γ_d1(machine)
    γ_q1 = PSY.get_γ_q1(machine)
    γ_d2 = PSY.get_γ_d2(machine)

    #Additional Fluxes
    ψd_pp = γ_d1 * eq_p + γ_q1 * ψ_kd
    ψ_pp = sqrt.(ψd_pp .^ 2 + ψq_pp .^ 2)

    I_fd = similar(δ, Float64)

    for ix in 1:length(δ)
        v = δ[ix]
        V_d, V_q = ri_dq(v) * [V_R[ix]; V_I[ix]]

        #Obtain electric current
        I_d =
            (1.0 / (R^2 + Xd_pp^2)) * (-R * (V_d - ψq_pp[ix]) + Xq_pp * (-V_q + ψd_pp[ix]))
        Se = saturation_function(machine, ψ_pp[ix])
        I_fd[ix] =
            eq_p[ix] +
            Se * ψd_pp[ix] +
            (Xd - Xd_p) * (I_d + γ_d2 * (eq_p[ix] - ψ_kd[ix] - (Xd_p - Xl) * I_d))
    end
    return ts, I_fd
end

"""
Function to obtain the field voltage time series of a Dynamic Generator with avrs that have 
the field voltage as a state. By default it is assumed that the models have that state.

"""
function _field_voltage(
    ::A,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
) where {A <: PSY.AVR}
    return post_proc_state_series(res, (name, :Vf), dt)
end

"""
Function to obtain the field voltage time series of a Dynamic Generator with avr AVRFixed.

"""
function _field_voltage(
    avr::PSY.AVRFixed,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    Vf0 = PSY.get_Vf(avr)
    ts, _ = post_proc_state_series(res, (name, :δ), dt)
    Vf = Vf0 * ones(length(ts))
    return ts, Vf
end

"""
Function to obtain the field voltage time series of a Dynamic Generator with avr ESAC1A and EXAC1.

"""
function _field_voltage(
    avr::Union{PSY.ESAC1A, PSY.EXAC1},
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, Ve = post_proc_state_series(res, (name, :Ve), dt)
    _, Xad_Ifd = post_proc_field_current_series(res, name, dt)
    Kc = PSY.get_Kc(avr)
    I_N = Kc * Xad_Ifd ./ Ve
    Vf = Ve .* rectifier_function.(I_N)
    return ts, Vf
end

"""
Function to obtain the field voltage time series of a Dynamic Generator with avr SCRX.

"""
function _field_voltage(
    avr::PSY.SCRX,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, Vr2 = post_proc_state_series(res, (name, :Vr2), dt)
    _, Ifd = post_proc_field_current_series(res, name, dt)
    V_min, V_max = PSY.get_Efd_lim(avr)
    bus_str = split(name, "-")[2]
    bus_num = parse(Int, bus_str)
    _, Vm = get_voltage_magnitude_series(res, bus_num; dt = dt)
    switch = PSY.get_switch(avr)
    rc_rfd = PSY.get_rc_rfd(avr)
    mult = switch == 0 ? Vm : ones(length(Vm))
    Vr2_sat = clamp.(Vr2, V_min, V_max)
    V_ex = Vr2_sat .* mult
    Vf = similar(V_ex)
    if rc_rfd == 0.0
        Vf = V_ex
    else
        for (ix, Ifd_ix) in enumerate(Ifd)
            if Ifd_ix > 0.0
                Vf[ix] = V_ex[ix]
            else
                Vf[ix] = -Ifd_ix * rc_rfd
            end
        end
    end
    return ts, Vf
end

"""
Function to obtain the field voltage time series of a Dynamic Generator with avr ESST1A.

"""
function _field_voltage(
    avr::PSY.ESST1A,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # Obtain state Va
    ts, Va = post_proc_state_series(res, (name, :Va), dt)

    # Obtain field current
    _, Ifd = post_proc_field_current_series(res, name, dt)

    # Obtain PSS output
    _, Vs = post_proc_pss_output_series(res, name, dt)

    # Get parameters
    VOS = PSY.get_PSS_flags(avr)
    Vr_min, Vr_max = PSY.get_Vr_lim(avr)
    Kc = PSY.get_Kc(avr)
    K_lr = PSY.get_K_lr(avr)
    I_lr = PSY.get_I_lr(avr)

    # Obtain machine's terminal voltage
    bus_str = split(name, "-")[2]
    bus_num = parse(Int, bus_str)
    _, Vt = get_voltage_magnitude_series(res, bus_num; dt = dt)

    # Compute auxiliary arrays
    Itemp = similar(Ifd)
    Iresult = similar(Ifd)

    for (ix, Ifd_ix) in enumerate(Ifd)
        Itemp[ix] = K_lr * (Ifd_ix - I_lr)
        Iresult[ix] = Itemp[ix] > 0.0 ? Itemp[ix] : 0.0
    end

    # Compute Va_sum
    Va_sum = similar(Va)

    if VOS == 1
        for (ix, Va_ix) in enumerate(Va)
            Va_sum[ix] = Va_ix - Iresult[ix]
        end
    elseif VOS == 2
        for (ix, Va_ix) in enumerate(Va)
            Va_sum[ix] = Va_ix - Iresult[ix] + Vs[ix]
        end
    end

    # Compute field voltage
    Vf = similar(Va_sum)
    for (ix, Va_sum_ix) in enumerate(Va_sum)
        Vf[ix] = clamp(Va_sum_ix, Vt[ix] * Vr_min, Vt[ix] * Vr_max - Kc * Ifd[ix])
    end

    return ts, Vf
end

"""
Function to obtain the pss output time series of a Dynamic Generator with pss PSSFixed.

"""
function _pss_output(
    pss::PSY.PSSFixed,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    ts, _ = post_proc_state_series(res, (name, :δ), dt)
    return ts, zeros(length(ts))
end

"""
Function to obtain the pss output time series of a Dynamic Generator with pss IEEEST.

"""
function _pss_output(
    pss::PSY.IEEEST,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # Obtain states
    ts, x_p2 = post_proc_state_series(res, (name, :x_p2), dt)
    _, x_p3 = post_proc_state_series(res, (name, :x_p3), dt)
    _, x_p4 = post_proc_state_series(res, (name, :x_p4), dt)
    _, x_p5 = post_proc_state_series(res, (name, :x_p5), dt)
    _, x_p6 = post_proc_state_series(res, (name, :x_p6), dt)
    _, x_p7 = post_proc_state_series(res, (name, :x_p7), dt)

    # Get Parameters
    A1 = PSY.get_A1(pss)
    A2 = PSY.get_A2(pss)
    A5 = PSY.get_A5(pss)
    A6 = PSY.get_A6(pss)
    T1 = PSY.get_T1(pss)
    T2 = PSY.get_T2(pss)
    T3 = PSY.get_T3(pss)
    T4 = PSY.get_T4(pss)
    T5 = PSY.get_T5(pss)
    T6 = PSY.get_T6(pss)
    Ks = PSY.get_Ks(pss)
    Ls_min, Ls_max = PSY.get_Ls_lim(pss)
    V_cu = PSY.get_Vcu(pss)
    V_cl = PSY.get_Vcl(pss)

    # Compute intermediate variables
    y_f = similar(x_p2)
    y_LL1 = similar(x_p2)
    y_LL2 = similar(x_p2)
    y_out = similar(x_p2)

    for (ix, _) in enumerate(x_p2)
        y_f[ix], _, _ =
            lead_lag_2nd_mass_matrix(x_p2[ix], x_p3[ix], x_p4[ix], A1, A2, A5, A6)
        y_LL1[ix], _ = lead_lag_mass_matrix(y_f[ix], x_p5[ix], 1.0, T1, T2)
        y_LL2[ix], _ = lead_lag_mass_matrix(y_LL1[ix], x_p6[ix], 1.0, T3, T4)
        y_out[ix], _ = high_pass_mass_matrix(y_LL2[ix], x_p7[ix], Ks * T5, T6)
    end

    # Compute V_ss
    V_ss = clamp.(y_out, Ls_min, Ls_max)

    # Obtain machine's terminal voltage
    bus_str = split(name, "-")[2]
    bus_num = parse(Int, bus_str)
    # TODO: Figure out how to compensate terminal voltage
    _, V_ct = get_voltage_magnitude_series(res, bus_num; dt = dt)

    # Compute PSS output signal
    V_pss = similar(V_ss)

    for (ix, V_ss_ix) in enumerate(V_ss)
        V_pss[ix] = output_pss_limiter(V_ss_ix, V_ct[ix], V_cl, V_cu)
    end

    return ts, V_pss
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with TGFixed Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.TGFixed,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # TODO: This will not plot correctly when changing P_ref in a callback
    ts, _ = _post_proc_state_series(res.solution, 1, dt)
    setpoints = get_setpoints(res)
    P_ref = setpoints[name]["P_ref"]
    efficiency = PSY.get_efficiency(tg)
    τm0 = P_ref * efficiency
    return ts, τm0 * ones(length(ts))
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with TGTypeI Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.TGTypeI,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # Get params
    Tc = PSY.get_Tc(tg)
    T3 = PSY.get_T3(tg)
    T4 = PSY.get_T4(tg)
    T5 = PSY.get_T5(tg)
    # Get state results
    ts, x_g1 = post_proc_state_series(res, (name, :x_g1), dt)
    _, x_g2 = post_proc_state_series(res, (name, :x_g2), dt)
    _, x_g3 = post_proc_state_series(res, (name, :x_g3), dt)
    τm = zeros(length(ts))
    for ix in 1:length(ts)
        y_ll, _ = lead_lag(x_g1[ix], x_g2[ix], 1.0, T3, Tc)
        τ_out, _ = lead_lag(y_ll, x_g3[ix], 1.0, T4, T5)
        τm[ix] = τ_out
    end
    return ts, τm
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with TGTypeII Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.TGTypeII,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # TODO: This will not plot correctly when changing P_ref in a callback
    # Get params
    setpoints = get_setpoints(res)
    P_ref = setpoints[name]["P_ref"]
    ω_ref = setpoints[name]["ω_ref"]
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    # Get state results
    ts, xg = post_proc_state_series(res, (name, :xg), dt)
    _, ω = post_proc_state_series(res, (name, :ω), dt)
    τm = zeros(length(ts))
    for ix in 1:length(ts)
        y_ll, _ = lead_lag(inv_R * (ω_ref - ω[ix]), xg[ix], 1.0, T1, T2)
        τm[ix] = y_ll + P_ref
    end
    return ts, τm
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with SteamTurbineGov1 (TGOV1) Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.SteamTurbineGov1,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # TODO: This will not plot correctly when changing P_ref in a callback
    # Get params
    setpoints = get_setpoints(res)
    P_ref = setpoints[name]["P_ref"]
    inv_R = PSY.get_R(tg) < eps() ? 0.0 : (1.0 / PSY.get_R(tg))
    T1 = PSY.get_T1(tg)
    T2 = PSY.get_T2(tg)
    V_min, V_max = PSY.get_valve_position_limits(tg)
    T3 = PSY.get_T3(tg)
    D_T = PSY.get_D_T(tg)
    # Get state results
    ts, x_g1 = post_proc_state_series(res, (name, :x_g1), dt)
    _, x_g2 = post_proc_state_series(res, (name, :x_g2), dt)
    _, ω = post_proc_state_series(res, (name, :ω), dt)
    ref_in = inv_R * (P_ref .- (ω .- 1.0))
    τm = zeros(length(ts))
    for ix in 1:length(ts)
        x_g1_sat, _ = low_pass_nonwindup(ref_in[ix], x_g1[ix], 1.0, T1, V_min, V_max)
        y_ll, _ = lead_lag(x_g1_sat, x_g2[ix], 1.0, T2, T3)
        P_m = y_ll - D_T * (ω[ix] - 1.0)
        τm[ix] = P_m / ω[ix]
    end
    return ts, τm
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with GasTG (GAST) Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.GasTG,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # Get params
    D_turb = PSY.get_D_turb(tg)
    # Get state results
    ts, x_g2 = post_proc_state_series(res, (name, :x_g2), dt)
    _, ω = post_proc_state_series(res, (name, :ω), dt)
    P_m = x_g2 - D_turb * (ω .- 1.0)
    τm = P_m ./ ω
    return ts, τm
end

"""
Function to obtain the mechanical torque time series of a Dynamic Generator with DEGOV Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.DEGOV,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64},
)
    ts, x_a3 = post_proc_state_series(res, (name, :x_a3), dt)
    return ts, x_a3
end
"""
Function to obtain the mechanical torque time series of a Dynamic Generator with HydroTurbineGov (HYGOV) Turbine Governor.

"""
function _mechanical_torque(
    tg::PSY.HydroTurbineGov,
    name::String,
    res::SimulationResults,
    dt::Union{Nothing, Float64, Vector{Float64}},
)
    # Get params
    q_nl = PSY.get_q_nl(tg)
    At = PSY.get_At(tg)
    D_T = PSY.get_D_T(tg)
    setpoints = get_setpoints(res)
    ω_ref = setpoints[name]["ω_ref"]
    # Get state results
    ts, x_g3 = post_proc_state_series(res, (name, :x_g3), dt)
    _, x_g4 = post_proc_state_series(res, (name, :x_g4), dt)
    _, ω = post_proc_state_series(res, (name, :ω), dt)
    Δω = ω .- ω_ref
    h = (x_g4 ./ x_g3) .^ 2
    τm = ((x_g4 .- q_nl) .* h * At - D_T * Δω .* x_g3) ./ ω
    return ts, τm
end
