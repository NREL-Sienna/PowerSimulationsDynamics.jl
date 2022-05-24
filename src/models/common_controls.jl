## Common implementation of each block:
## returns a tuple of (output, internal_derivative)

# Low pass modified with denominator K_den instead of 1.0
"""
Low Pass Filter Modified
     ┌─────────────┐
     │      K      │
u -> │ ────────────│ -> y
     │ K_den + sT  │
     └─────────────┘
"""
function low_pass_modified_mass_matrix(
    u::Z,
    y::V,
    K::Float64,
    K_den::W,
    ::Float64,
) where {V <: ACCEPTED_REAL_TYPES, W <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    return y, K * u - K_den * y
end

function low_pass_modified(
    u::Z,
    y::V,
    K::Float64,
    K_den::W,
    T::Float64,
) where {V <: ACCEPTED_REAL_TYPES, W <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    return y, (1.0 / T) * low_pass_modified_mass_matrix(u, y, K, K_den, T)[2]
end

"""
Low Pass Filter
     ┌────────┐
     │    K   │
u -> │ ────── │ -> y
     │ 1 + sT │
     └────────┘
"""

# Use this one if T = 0 is allowed, and let the mass matrix take care of it.
function low_pass_mass_matrix(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    return low_pass_modified_mass_matrix(u, y, K, 1.0, T)
end

function low_pass(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    return low_pass_modified(u, y, K, 1.0, T)
end

"""
Low Pass Filter with non-windup
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

function low_pass_nonwindup_mass_matrix(
    u::Z,
    y::V,
    K::Float64,
    ::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    dydt = K * u - y
    y_sat = clamp(y, y_min, y_max)
    # Non Windup logic from IEEE Std 421.5
    binary_logic = ((y >= y_max) && (dydt > 0)) || ((y <= y_min) && (dydt < 0)) ? 0.0 : 1.0
    return y_sat, binary_logic * dydt
end

# Does not accept T = 0
function low_pass_nonwindup(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    y_sat, dydt_scaled = low_pass_nonwindup_mass_matrix(u, y, K, T, y_min, y_max)
    return y_sat, (1.0 / T) * dydt_scaled
end

# Does not accept T = 0
function low_pass_nonwindup_ramp_limits(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
    dy_min::Float64,
    dy_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    y, dydt = low_pass(u, y, K, T)
    y_sat = clamp(y, y_min, y_max)
    dydt_sat = clamp(dydt, dy_min, dy_max)
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

function high_pass_mass_matrix(
    u::Z,
    x::Z,
    K::Float64,
    T::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    K_T = T < eps() ? 0.0 : (K / T)
    return x + K_T * u, -(K_T * u + x)
end

function high_pass(u::Z, x::Z, K::Float64, T::Float64) where {Z <: ACCEPTED_REAL_TYPES}
    y, dxdt_scaled = high_pass_mass_matrix(u, x, K, T)
    return y, (1.0 / T) * dxdt_scaled
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

function lead_lag_mass_matrix(
    u::Z,
    x::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    T1_T2 = T2 < eps() ? 0.0 : (T1 / T2)
    return x + (K * T1_T2) * u, K * (1 - T1_T2) * u - x
end

function lead_lag(
    u::Z,
    x::Z,
    K::Float64,
    T1::Float64,
    T2::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dxdt_scaled = lead_lag_mass_matrix(u, x, K, T1, T2)
    return y, (1.0 / T2) * dxdt_scaled
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
    dxdt = K * u - T1 * x - y
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

function lead_lag_2nd_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T1 * x1 - x2
    dx2dt = x1
    T4_T2 = T2 < eps() ? 0.0 : (T4 / T2)
    y = T4_T2 * u + (T3 - T1 * T4_T2) * x1 + (1.0 - T4_T2) * x2
    return y, dx1dt, dx2dt
end

function lead_lag_2nd(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt = lead_lag_2nd_mass_matrix(u, x1, x2, T1, T2, T3, T4)
    return y, (1.0 / T2) * dx1dt_scaled, dx2dt
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
    y, _ = pi_block(u, x, kp, ki)
    y_sat = clamp(y, y_min, y_max)
    binary_logic = y_min < y < y_max ? 1.0 : 0.0
    return y_sat, binary_logic * u
end

"""
Integrator with windup limits
                                 y_max
                                _ _ _ 
        ┌────────┐             /
        │    K   │   y        /
 u - -->│ ────── │ - - - - - /- - - - - --> ysat
        │   sT   │          /
        └────────┘   _ _ _ / 
                     y_min
 """

function integrator_windup_mass_matrix(
    u::Z,
    y::V,
    K::Float64,
    ::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    dydt_scaled = K * u
    y_sat = clamp(y, y_min, y_max)
    return y_sat, dydt_scaled
end

# Does not accept T = 0
function integrator_windup(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    y_sat = clamp(y, y_min, y_max)
    return y_sat, (1.0 / T) * integrator_windup_mass_matrix(u, y, K, T, y_min, y_max)[2]
end

"""
Integrator with non-windup limits
             y_max
           /¯¯¯¯¯¯
     ┌────────┐
     │    K   │
u -> │ ────── │ -> y
     │   sT   │
     └────────┘
   ______/
   y_min
"""

function integrator_nonwindup_mass_matrix(
    u::Z,
    y::V,
    K::Float64,
    ::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    dydt_scaled = K * u
    y_sat = clamp(y, y_min, y_max)
    binary_logic =
        ((y >= y_max) && (dydt_scaled > 0)) || ((y <= y_min) && (dydt_scaled < 0)) ? 0.0 :
        1.0
    return y_sat, binary_logic * dydt_scaled
end

# Does not accept T = 0
function integrator_nonwindup(
    u::Z,
    y::V,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {V <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    y_sat, dydt_scaled = integrator_nonwindup_mass_matrix(u, y, K, T, y_min, y_max)
    return y_sat, (1.0 / T) * dydt_scaled
end
