"""
    Saturation function for quadratic saturation models for machines
        Se(x) = B * (x - A)^2 / x
"""
function saturation_function(
    machine::Union{PSY.RoundRotorQuadratic, PSY.SalientPoleQuadratic},
    x::Number,
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
    x::Number,
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

function saturation_function(avr::PSY.ESAC1A, x::Number)
    Sat_A, Sat_B = PSY.get_saturation_coeffs(avr)
    return Sat_B * (x - Sat_A)^2 / x
end

function rectifier_function(I::T) where {T <: Number}
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
