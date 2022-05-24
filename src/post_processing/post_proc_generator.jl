"""
Function to obtain the output current time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
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
Function to obtain the field current time series of a Dynamic Generator model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_field_current(
    res::SimulationResults,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
) where {G <: PSY.DynamicGenerator}

    #Get AVR
    avr = PSY.get_avr(dynamic_device)
    return _field_voltage(avr, PSY.get_name(dynamic_device), res, dt)
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
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
    dt::Union{Nothing, Float64},
)
    ts, Ve = post_proc_state_series(res, (name, :Ve), dt)
    _, Xad_Ifd = post_proc_field_current_series(res, name, dt)
    Kc = PSY.get_Kc(avr)
    I_N = Kc * Xad_Ifd ./ Ve
    Vf = Ve .* rectifier_function.(I_N)
    return ts, Vf
end
