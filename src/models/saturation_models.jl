"""
    Saturation function for quadratic saturation models for machines
        Se(x) = B * (x - A)^2 / x
"""
function saturation_function(
    machine::Union{PSY.RoundRotorQuadratic, PSY.SalientPoleQuadratic},
    x::Real,
)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(machine)
    return Sat_B * (x - Sat_A)^2 / x
end

"""
    Saturation function for exponential saturation models for machines
        Se(x) = B * x^A
"""
function saturation_function(
    machine::Union{PSY.RoundRotorExponential, PSY.SalientPoleExponential},
    x::Real,
)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(machine)
    return Sat_B * x^Sat_A
end

function rectifier_function(I::Float64)
    if I <= 0.0
        return 1.0
    elseif I <= 0.433
        return 1.0 - 0.577 * I
    elseif I < 0.75
        return sqrt(0.75 - I^2)
    elseif I <= 1.0
        return 1.732 * (1.0 - I)
    else
        return 0.0
    end
end

function saturation_function(avr::PSY.ESAC1A, x::Real)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(avr)
    return Sat_B * (x - Sat_A)^2 / x
end

function rectifier_function(I::T) where {T <: Real}
    if I <= 0.0
        return one(T)
    elseif I <= 0.433
        return 1.0 - 0.577 * I
    elseif I < 0.75
        return sqrt(0.75 - I^2)
    elseif I <= 1.0
        return 1.732 * (1.0 - I)
    else
        return zero(T)
    end
end

function deadband_function(x::T, db_low::Float64, db_high::Float64) where {T <: Real}
    if x > db_high
        return x - db_high
    elseif x < db_low
        return x - db_low
    else
        return zero(T)
    end
end

function current_limit_logic(inner_control::InnerREECB1, Vt_filt::T, Ip_cmd::T, Iq_cmd::T) where {T <: Real}
    I_max = PSY.get_I_max(inner_control)
    PQ_Flag = PSY.get_PQ_Flag(inner_control)
    Iq_max = I_max
    Ip_max = I_max

    if PQ_Flag == 0 # Q Priority
        Iq_min = - Iq_max
        local_I = I_max^2 - Iq_cmd^2
        if local_I < 0
            local_I = 0
        else
            local_I = sqrt(local_I)
        end
        if local_I < Ip_max
            Ip_max = local_I
        end
        Ip_min = 0
    else # P Priority
        Ip_min = 0
        local_I = I_max^2 - Ip_cmd^2
        if local_I < 0
            local_I = 0
        else
            local_I = sqrt(local_I)
        end
        if local_I < Iq_max
            Iq_max = local_I
        end
        Iq_min = - Iq_max
    end
    return Ip_min, Ip_max, Iq_min, Iq_max
end
