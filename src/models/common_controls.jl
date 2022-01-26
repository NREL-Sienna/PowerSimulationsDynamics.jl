## Common implementation of each block:
## returns a tuple of (output, internal_derivative)

"""
Low Pass Filter
             y_max
           /¯¯¯¯¯¯
     ┌────────┐
     │    K   │
u -> │ ────── │ -> y
     │ 1 + sT │
     └────────┘
   ______/
   y_min
"""

function low_pass(u::Z, y::Z, K::Float64, T::Float64) where {Z <: ACCEPTED_REAL_TYPES}
    return y, (1.0 / T) * (K * u - y)
end

# Low pass modified with denominator K_den instead of 1.0
"""
Low Pass Filter Modified
     ┌─────────────┐
     │      K      │
u -> │ ────────────│ -> y
     │ K_den + sT  │
     └─────────────┘
"""
function low_pass_modified(
    u::Z,
    y::Z,
    K::Float64,
    K_den::Z,
    T::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    return y, (1.0 / T) * (K * u - K_den * y)
end

# Use this one if T = 0 is allowed, and let the mass matrix take care of it.
function low_pass_mass_matrix(
    u::Z,
    y::Z,
    K::Float64,
    T::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    return y, T * low_pass(u, y, K, T)[2]
end

# Does not accept T = 0
function low_pass_nonwindup(
    u::Z,
    y::Z,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dydt = (1.0 / T) * (K * u - y)
    y_sat = clamp(y, y_min, y_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic = ((y >= y_max) && (dydt > 0)) || ((y <= y_min) && (dydt < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt
end

function low_pass_nonwindup_mass_matrix(
    u::Z,
    y::Z,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dydt = K * u - y
    y_sat = clamp(y, y_min, y_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic = ((y >= y_max) && (dydt > 0)) || ((y <= y_min) && (dydt < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt
end

# Does not accept T = 0
function low_pass_nonwindup_ramp_limits(
    u::Z,
    y::Z,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
    dy_min::Float64,
    dy_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y_sat = clamp(y, y_min, y_max)
    dydt_sat = clamp((1.0 / T) * (K * u - y), dy_min, dy_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic =
        ((y >= y_max) && (dydt_sat > 0)) || ((y <= y_min) && (dydt_sat < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt_sat
end

"""
High Pass Filter 
     ┌────────┐
     │   sK   │
u -> │ ────── │ -> y
     │ 1 + sT │
     └────────┘

Internal State: x
"""

function high_pass(u::Z, x::Z, K::Float64, T::Float64) where {Z <: ACCEPTED_REAL_TYPES}
    dxdt = -(1.0 / T) * ((K / T) * u + x)
    return x + (K / T) * u, dxdt
end

function high_pass_mass_matrix(
    u::Z,
    x::Z,
    K::Float64,
    T::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    return x + (K / T) * u, T * high_pass(u, x, K, T)[2]
end

"""
Lead-Lag Block:
     ┌───────────┐
     │   1 + sT1 │
u -> │ K ─────── │ -> y
     │   1 + sT2 │
     └───────────┘

Internal State: x
"""
function lead_lag(
    u::Z,
    x::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dxdt = (1.0 / T2) * (K * (1 - T1 / T2) * u - x)
    return x + (K * T1 / T2) * u, dxdt
end

function lead_lag_mass_matrix(
    u::Z,
    x::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    return x + (K * T1 / T2) * u, T2 * lead_lag(u, x, K, T1, T2)[2]
end

"""
2nd Order Low Pass Filter
     ┌──────────────────┐
     │         K        │
u -> │ ──────────────── │ -> y
     │ 1 + sT1 + s^2 T2 │
     └──────────────────┘

Internal state: x
"""

function low_pass_2nd(
    u::Z,
    x::Z,
    y::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dxdt = (1.0 / T2) * (K * u - T1 * x - y)
    dydt = x
    return y, dxdt, dydt
end

function low_pass_2nd_mass_matrix(
    u::Z,
    x::Z,
    y::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dxdt = T2 * low_pass_2nd(u, x, y, K, T1, T2)[2]
    dydt = x
    return y, dxdt, dydt
end

"""
2nd Order Lead-Lag Block
     ┌──────────────────┐
     │ 1 + sT3 + s^2 T4 │
u -> │ ──────────────── │ -> y
     │ 1 + sT1 + s^2 T2 │
     └──────────────────┘
    
Internal States: x1, x2
"""

function lead_lag_2nd(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = x2
    dx2dt = (1.0 / T2) * (u - T1 * x2 - x1)
    y = (T4 / T2) * u + (T3 - T4 * T1 / T2) * x2 + (1.0 - T4 / T2) * x1
    return y, dx1dt, dx2dt
end

function lead_lag_2nd_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = x2
    dx2dt = u - T1 * x2 - x1
    y = (T4 / T2) * u + (T3 - T4 * T1 / T2) * x2 + (1.0 - T4 / T2) * x1
    return y, dx1dt, dx2dt
end

"""
Proportional-Integral Block
             y_max
            /¯¯¯¯¯¯
     ┌──────────┐
     │      ki  │
u -> │kp + ───  │ -> y
     │      s   │
     └──────────┘
    ______/
     y_min

Internal State: x
"""

function pi_block(u::Z, x::Z, kp::Float64, ki::Float64) where {Z <: ACCEPTED_REAL_TYPES}
    return kp * u + ki * x, u
end

function pi_block_nonwindup(
    u::Z,
    x::Z,
    kp::Float64,
    ki::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y = kp * u + ki * x
    y_sat = clamp(y, y_min, y_max)
    binary_logic = y_min < y < y_max ? 1.0 : 0.0
    return y_sat, binary_logic * u
end
