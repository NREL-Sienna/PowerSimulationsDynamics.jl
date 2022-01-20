## Common implementation of each block:
## returns a tuple of (output, internal_derivative)

function low_pass(
    u::X,
    y::X,
    K::Float64,
    T::Float64
) where {X <: ACCEPTED_REAL_TYPES}
    return y, (1.0 / T) * (K * u - y)
end

# Use this one if T = 0 is allowed, and let the mass matrix take care of it.
function low_pass_mass_matrix(
    u::X,
    y::X,
    K::Float64,
    T::Float64
) where {X <: ACCEPTED_REAL_TYPES}
    return y, T * low_pass(u, y, K, T)[2]
end

# Does not accept T = 0
function low_pass_nonwindup(
    u::X,
    y::X,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {X <: ACCEPTED_REAL_TYPES}
    dydt = (1.0 / T) * (K * u - y)
    y_sat = clamp(y, y_min, y_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic = ((y >= y_max) && (dydt > 0)) || ((y <= y_min) && (dydt < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt
end

# Does not accept T = 0
function low_pass_nonwindup_ramp_limits(
    u::X,
    y::X,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
    dy_min::Float64,
    dy_max::Float64,
) where {X <: ACCEPTED_REAL_TYPES}
    y_sat = clamp(y, y_min, y_max)
    dydt_sat = clamp((1.0 / T) * (K * u - y), dy_min, dy_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic = ((y >= y_max) && (dydt_sat > 0)) || ((y <= y_min) && (dydt_sat < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt_sat
end