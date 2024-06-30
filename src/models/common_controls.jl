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
    K::ACCEPTED_REAL_TYPES,
    K_den::W,
    ::ACCEPTED_REAL_TYPES,
) where {V <: ACCEPTED_REAL_TYPES, W <: ACCEPTED_REAL_TYPES, Z <: ACCEPTED_REAL_TYPES}
    return y, K * u - K_den * y
end

function low_pass_modified(
    u::Z,
    y::V,
    K::ACCEPTED_REAL_TYPES,
    K_den::W,
    T::ACCEPTED_REAL_TYPES,
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
    K::ACCEPTED_REAL_TYPES,
    T::ACCEPTED_REAL_TYPES,
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
    K::ACCEPTED_REAL_TYPES,
    T::ACCEPTED_REAL_TYPES,
) where {Z <: ACCEPTED_REAL_TYPES}
    K_T = T < eps() ? 0.0 : (K / T)
    return x + K_T * u, -(K_T * u + x)
end

function high_pass(
    u::Z,
    x::Z,
    K::ACCEPTED_REAL_TYPES,
    T::ACCEPTED_REAL_TYPES,
) where {Z <: ACCEPTED_REAL_TYPES}
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
    T1::ACCEPTED_REAL_TYPES,
    T2::ACCEPTED_REAL_TYPES,
    T3::ACCEPTED_REAL_TYPES,
    T4::ACCEPTED_REAL_TYPES,
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
2nd Order Lead-Lag Block with non-windup limits
                  y_max
                 /¯¯¯¯¯¯
     ┌──────────────────┐
     │ 1 + sT3 + s^2 T4 │
u -> │ ──────────────── │ -> y
     │ 1 + sT1 + s^2 T2 │
     └──────────────────┘
      ______/
      y_min
    
Internal States: x1, x2
"""

function lead_lag_2nd_mass_matrix_nonwindup(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt_scaled = u - T1 * x1 - x2
    dx2dt = x1
    T4_T2 = T2 < eps() ? 0.0 : (T4 / T2)
    y = T4_T2 * u + (T3 - T1 * T4_T2) * x1 + (1.0 - T4_T2) * x2
    y_sat = clamp(y, y_min, y_max)
    binary_logic = y_min < y < y_max ? 1.0 : 0.0
    return y_sat, binary_logic * dx1dt_scaled, binary_logic * dx2dt
end

function lead_lag_2nd_nonwindup(
    u::Z,
    x1::Z,
    x2::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt =
        lead_lag_2nd_mass_matrix_nonwindup(u, x1, x2, T1, T2, T3, T4, y_min, y_max)
    return y, (1.0 / T2) * dx1dt_scaled, dx2dt
end

"""
3rd Order Lead-Lag Block
     ┌───────────────────────────┐
     │ 1 + sT4 + s^2 T5 + s^3 T6 │
u -> │ ───────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2 + s^3 T3 │
     └───────────────────────────┘
    
Internal States: x1, x2, x3
"""

function lead_lag_3rd_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T2 * x1 - T1 * x2 - x3
    dx2dt = x1
    dx3dt = x2
    T6_T3 = T3 < eps() ? 0.0 : (T6 / T3)
    y = T6_T3 * u + (T5 - T2 * T6_T3) * x1 + (T4 - T1 * T6_T3) * x2 + (1.0 - T6_T3) * x3
    return y, dx1dt, dx2dt, dx3dt
end

function lead_lag_3rd(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt =
        lead_lag_3rd_mass_matrix(u, x1, x2, x3, T1, T2, T3, T4, T5, T6)
    return y, (1.0 / T3) * dx1dt_scaled, dx2dt, dx3dt
end

"""
4th Order Lead-Lag Block
     ┌────────────────────────────────────┐
     │ 1 + sT5 + s^2 T6 + s^3 T7 + s^4 T8 │
u -> │ ────────────────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2 + s^3 T3 + s^4 T4 │
     └────────────────────────────────────┘
    
Internal States: x1, x2, x3, x4
"""

function lead_lag_4th_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T3 * x1 - T2 * x2 - T1 * x3 - x4
    dx2dt = x1
    dx3dt = x2
    dx4dt = x3
    T8_T4 = T4 < eps() ? 0.0 : (T8 / T4)
    y =
        T8_T4 * u +
        (T7 - T3 * T8_T4) * x1 +
        (T6 - T2 * T8_T4) * x2 +
        (T5 - T1 * T8_T4) * x3 +
        (1.0 - T8_T4) * x4
    return y, dx1dt, dx2dt, dx3dt, dx4dt
end

function lead_lag_4th(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt, dx4dt =
        lead_lag_4th_mass_matrix(u, x1, x2, x3, x4, T1, T2, T3, T4, T5, T6, T7, T8)
    return y, (1.0 / T4) * dx1dt_scaled, dx2dt, dx3dt, dx4dt
end

"""
5th Order Lead-Lag Block
     ┌──────────────────────────────────────────────┐
     │ 1 + sT6 + s^2 T7 + s^3 T8 + s^4 T9 + s^5 T10 │
u -> │ ──────────────────────────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2 + s^3 T3 + s^4 T4 + s^5 T5  │
     └──────────────────────────────────────────────┘
    
Internal States: x1, x2, x3, x4, x5
"""

function lead_lag_5th_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T4 * x1 - T3 * x2 - T2 * x3 - T1 * x4 - x5
    dx2dt = x1
    dx3dt = x2
    dx4dt = x3
    dx5dt = x4
    T10_T5 = T5 < eps() ? 0.0 : (T10 / T5)
    y =
        T10_T5 * u +
        (T9 - T4 * T10_T5) * x1 +
        (T8 - T3 * T10_T5) * x2 +
        (T7 - T2 * T10_T5) * x3 +
        (T6 - T1 * T10_T5) * x4 +
        (1.0 - T10_T5) * x5
    return y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt
end

function lead_lag_5th(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt = lead_lag_5th_mass_matrix(
        u,
        x1,
        x2,
        x3,
        x4,
        x5,
        T1,
        T2,
        T3,
        T4,
        T5,
        T6,
        T7,
        T8,
        T9,
        T10,
    )
    return y, (1.0 / T5) * dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt
end

"""
6th Order Lead-Lag Block
     ┌─────────────────────────────────────────────────────────┐
     │ 1 + sT7 + s^2 T8 + s^3 T9 + s^4 T10 + s^5 T11 + s^6 T12 │
u -> │ ─────────────────────────────────────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2 + s^3 T3 + s^4 T4  + s^5 T5  + s^6 T6  │
     └─────────────────────────────────────────────────────────┘
    
Internal States: x1, x2, x3, x4, x5, x6
"""

function lead_lag_6th_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T5 * x1 - T4 * x2 - T3 * x3 - T2 * x4 - T1 * x5 - x6
    dx2dt = x1
    dx3dt = x2
    dx4dt = x3
    dx5dt = x4
    dx6dt = x5
    T12_T6 = T6 < eps() ? 0.0 : (T12 / T6)
    y =
        T12_T6 * u +
        (T11 - T5 * T12_T6) * x1 +
        (T10 - T4 * T12_T6) * x2 +
        (T9 - T3 * T12_T6) * x3 +
        (T8 - T2 * T12_T6) * x4 +
        (T7 - T1 * T12_T6) * x5 +
        (1.0 - T12_T6) * x6
    return y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt
end

function lead_lag_6th(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt = lead_lag_6th_mass_matrix(
        u,
        x1,
        x2,
        x3,
        x4,
        x5,
        x6,
        T1,
        T2,
        T3,
        T4,
        T5,
        T6,
        T7,
        T8,
        T9,
        T10,
        T11,
        T12,
    )
    return y, (1.0 / T6) * dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt
end

"""
7th Order Lead-Lag Block
     ┌────────────────────────────────────────────────────────────────────┐
     │ 1 + sT8 + s^2 T9 + s^3 T10 + s^4 T11 + s^5 T12 + s^6 T13 + s^7 T14 │
u -> │ ────────────────────────────────────────────────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2 + s^3 T3  + s^4 T4  + s^5 T5  + s^6 T6  + s^7 T7  │
     └────────────────────────────────────────────────────────────────────┘
    
Internal States: x1, x2, x3, x4, x5, x6, x7
"""

function lead_lag_7th_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    x7::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
    T13::Float64,
    T14::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T6 * x1 - T5 * x2 - T4 * x3 - T3 * x4 - T2 * x5 - T1 * x6 - x7
    dx2dt = x1
    dx3dt = x2
    dx4dt = x3
    dx5dt = x4
    dx6dt = x5
    dx7dt = x6
    T14_T7 = T7 < eps() ? 0.0 : (T14 / T7)
    y =
        T14_T7 * u +
        (T13 - T6 * T14_T7) * x1 +
        (T12 - T5 * T14_T7) * x2 +
        (T11 - T4 * T14_T7) * x3 +
        (T10 - T3 * T14_T7) * x4 +
        (T9 - T2 * T14_T7) * x5 +
        (T8 - T1 * T14_T7) * x6 +
        (1.0 - T14_T7) * x7
    return y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt
end

function lead_lag_7th(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    x7::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
    T13::Float64,
    T14::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt = lead_lag_7th_mass_matrix(
        u,
        x1,
        x2,
        x3,
        x4,
        x5,
        x6,
        x7,
        T1,
        T2,
        T3,
        T4,
        T5,
        T6,
        T7,
        T8,
        T9,
        T10,
        T11,
        T12,
        T13,
        T14,
    )
    return y, (1.0 / T7) * dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt
end

"""
8th Order Lead-Lag Block
     ┌───────────────────────────────────────────────────────────────────────────────┐
     │ 1 + sT9 + s^2 T10 + s^3 T11 + s^4 T12 + s^5 T13 + s^6 T14 + s^7 T15 + s^8 T16 │
u -> │ ───────────────────────────────────────────────────────────────────────────── │ -> y
     │ 1 + sT1 + s^2 T2  + s^3 T3  + s^4 T4  + s^5 T5  + s^6 T6  + s^7 T7  + s^8 T8  │
     └───────────────────────────────────────────────────────────────────────────────┘
    
Internal States: x1, x2, x3, x4, x5, x6, x7, x8
"""

function lead_lag_8th_mass_matrix(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    x7::Z,
    x8::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
    T13::Float64,
    T14::Float64,
    T15::Float64,
    T16::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = u - T7 * x1 - T6 * x2 - T5 * x3 - T4 * x4 - T3 * x5 - T2 * x6 - T1 * x7 - x8
    dx2dt = x1
    dx3dt = x2
    dx4dt = x3
    dx5dt = x4
    dx6dt = x5
    dx7dt = x6
    dx8dt = x7
    T16_T8 = T8 < eps() ? 0.0 : (T16 / T8)
    y =
        T16_T8 * u +
        (T15 - T7 * T16_T8) * x1 +
        (T14 - T6 * T16_T8) * x2 +
        (T13 - T5 * T16_T8) * x3 +
        (T12 - T4 * T16_T8) * x4 +
        (T11 - T3 * T16_T8) * x5 +
        (T10 - T2 * T16_T8) * x6 +
        (T9 - T1 * T16_T8) * x7 +
        (1.0 - T16_T8) * x8
    return y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt
end

function lead_lag_8th(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    x7::Z,
    x8::Z,
    T1::Float64,
    T2::Float64,
    T3::Float64,
    T4::Float64,
    T5::Float64,
    T6::Float64,
    T7::Float64,
    T8::Float64,
    T9::Float64,
    T10::Float64,
    T11::Float64,
    T12::Float64,
    T13::Float64,
    T14::Float64,
    T15::Float64,
    T16::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y, dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt =
        lead_lag_8th_mass_matrix(
            u,
            x1,
            x2,
            x3,
            x4,
            x5,
            x6,
            x7,
            x8,
            T1,
            T2,
            T3,
            T4,
            T5,
            T6,
            T7,
            T8,
            T9,
            T10,
            T11,
            T12,
            T13,
            T14,
            T15,
            T16,
        )
    return y, (1.0 / T8) * dx1dt_scaled, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt
end

"""
Ramp Tracking Filter Block

M >= 0

N >= 0

M * N <= 8

If M == 0, N = 0

To bypass: use M = 0, N = 0

     ┌─────────────────────┐
     │  ┌─           ─┐^N  │
     │  │   1 + sT2   │    │
u -> │  │ ─────────── │    │ -> y
     │  │ (1 + sT1)^M │    │
     │  └─           ─┘    │
     └─────────────────────┘
    
Internal States: x1, x2, x3, x4, x5, x6, x7, x8
"""

function ramp_tracking_filter(
    u::Z,
    x1::Z,
    x2::Z,
    x3::Z,
    x4::Z,
    x5::Z,
    x6::Z,
    x7::Z,
    x8::Z,
    T1::Float64,
    T2::Float64,
    M::Int64,
    N::Int64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dx1dt = 0.0
    dx2dt = 0.0
    dx3dt = 0.0
    dx4dt = 0.0
    dx5dt = 0.0
    dx6dt = 0.0
    dx7dt = 0.0
    dx8dt = 0.0
    y = 0.0

    if N == 0
        y = u
    elseif N == 1
        if M == 0
            y = u
        elseif M == 1
            T1_ll = T2 # time parameter for the lead-lag filter
            T2_ll = T1 # time parameter for the lead-lag filter

            y, dx1dt = lead_lag(u, x1, 1.0, T1_ll, T2_ll)
        elseif M == 2
            T1_ll = 2.0 * T1 # time parameter for the lead-lag filter
            T2_ll = T1^2 # time parameter for the lead-lag filter
            T3_ll = T2 # time parameter for the lead-lag filter
            T4_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt = lead_lag_2nd(u, x1, x2, T1_ll, T2_ll, T3_ll, T4_ll)
        elseif M == 3
            T1_ll = 3.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 3.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = T1^3 # time parameter for the lead-lag filter
            T4_ll = T2 # time parameter for the lead-lag filter
            T5_ll = 0.0 # time parameter for the lead-lag filter
            T6_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt =
                lead_lag_3rd(u, x1, x2, x3, T1_ll, T2_ll, T3_ll, T4_ll, T5_ll, T6_ll)
        elseif M == 4
            T1_ll = 4.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 6.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 4.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = T1^4 # time parameter for the lead-lag filter
            T5_ll = T2 # time parameter for the lead-lag filter
            T6_ll = 0.0 # time parameter for the lead-lag filter
            T7_ll = 0.0 # time parameter for the lead-lag filter
            T8_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt = lead_lag_4th(
                u,
                x1,
                x2,
                x3,
                x4,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
            )
        elseif M == 5
            T1_ll = 5.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 10.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 10.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 5.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = T1^5 # time parameter for the lead-lag filter
            T6_ll = T2 # time parameter for the lead-lag filter
            T7_ll = 0.0 # time parameter for the lead-lag filter
            T8_ll = 0.0 # time parameter for the lead-lag filter
            T9_ll = 0.0 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt = lead_lag_5th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
            )
        elseif M == 6
            T1_ll = 6.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 15.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 20.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 15.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 6.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = T1^6 # time parameter for the lead-lag filter
            T7_ll = T2 # time parameter for the lead-lag filter
            T8_ll = 0.0 # time parameter for the lead-lag filter
            T9_ll = 0.0 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt = lead_lag_6th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
            )
        elseif M == 7
            T1_ll = 7.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 21.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 35.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 35.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 21.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 7.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = T1^7 # time parameter for the lead-lag filter
            T8_ll = T2 # time parameter for the lead-lag filter
            T9_ll = 0.0 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter
            T13_ll = 0.0 # time parameter for the lead-lag filter
            T14_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt = lead_lag_7th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
            )
        elseif M == 8
            T1_ll = 8.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 28.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 56.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 70.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 56.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 28.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = 8.0 * T1^7 # time parameter for the lead-lag filter
            T8_ll = T1^8 # time parameter for the lead-lag filter
            T9_ll = T2 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter
            T13_ll = 0.0 # time parameter for the lead-lag filter
            T14_ll = 0.0 # time parameter for the lead-lag filter
            T15_ll = 0.0 # time parameter for the lead-lag filter
            T16_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt = lead_lag_8th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                x8,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
                T15_ll,
                T16_ll,
            )
        end
    elseif N == 2
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 2.0 * T1 # time parameter for the lead-lag filter
            T2_ll = T1^2 # time parameter for the lead-lag filter
            T3_ll = 2.0 * T2 # time parameter for the lead-lag filter
            T4_ll = T2^2 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt = lead_lag_2nd(u, x1, x2, T1_ll, T2_ll, T3_ll, T4_ll)
        elseif M == 2
            T1_ll = 4.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 6.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 4.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = T1^4 # time parameter for the lead-lag filter
            T5_ll = 2.0 * T2 # time parameter for the lead-lag filter
            T6_ll = T2^2 # time parameter for the lead-lag filter
            T7_ll = 0.0 # time parameter for the lead-lag filter
            T8_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt = lead_lag_4th(
                u,
                x1,
                x2,
                x3,
                x4,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
            )
        elseif M == 3
            T1_ll = 6.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 15.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 20.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 15.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 6.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = T1^6 # time parameter for the lead-lag filter
            T7_ll = 2.0 * T2 # time parameter for the lead-lag filter
            T8_ll = T2^2 # time parameter for the lead-lag filter
            T9_ll = 0.0 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt = lead_lag_6th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
            )
        elseif M == 4
            T1_ll = 8.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 28.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 56.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 70.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 56.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 28.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = 8.0 * T1^7 # time parameter for the lead-lag filter
            T8_ll = T1^8 # time parameter for the lead-lag filter
            T9_ll = 2.0 * T2 # time parameter for the lead-lag filter
            T10_ll = T2^2 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter
            T13_ll = 0.0 # time parameter for the lead-lag filter
            T14_ll = 0.0 # time parameter for the lead-lag filter
            T15_ll = 0.0 # time parameter for the lead-lag filter
            T16_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt = lead_lag_8th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                x8,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
                T15_ll,
                T16_ll,
            )
        end
    elseif N == 3
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 3.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 3.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = T1^3 # time parameter for the lead-lag filter
            T4_ll = 3.0 * T2 # time parameter for the lead-lag filter
            T5_ll = 3.0 * T2^2 # time parameter for the lead-lag filter
            T6_ll = T2^3 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt =
                lead_lag_3rd(u, x1, x2, x3, T1_ll, T2_ll, T3_ll, T4_ll, T5_ll, T6_ll)
        elseif M == 2
            T1_ll = 6.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 15.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 20.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 15.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 6.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = T1^6 # time parameter for the lead-lag filter
            T7_ll = 3.0 * T2 # time parameter for the lead-lag filter
            T8_ll = 3.0 * T2^2 # time parameter for the lead-lag filter
            T9_ll = T2^3 # time parameter for the lead-lag filter
            T10_ll = 0.0 # time parameter for the lead-lag filter
            T11_ll = 0.0 # time parameter for the lead-lag filter
            T12_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt = lead_lag_6th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
            )
        end
    elseif N == 4
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 4.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 6.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 4.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = T1^4 # time parameter for the lead-lag filter
            T5_ll = 4.0 * T2 # time parameter for the lead-lag filter
            T6_ll = 6.0 * T2^2 # time parameter for the lead-lag filter
            T7_ll = 4.0 * T2^3 # time parameter for the lead-lag filter
            T8_ll = T2^4 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt = lead_lag_4th(
                u,
                x1,
                x2,
                x3,
                x4,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
            )
        elseif M == 2
            T1_ll = 8.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 28.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 56.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 70.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 56.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 28.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = 8.0 * T1^7 # time parameter for the lead-lag filter
            T8_ll = T1^8 # time parameter for the lead-lag filter
            T9_ll = 4.0 * T2 # time parameter for the lead-lag filter
            T10_ll = 6.0 * T2^2 # time parameter for the lead-lag filter
            T11_ll = 4.0 * T2^3 # time parameter for the lead-lag filter
            T12_ll = T2^4 # time parameter for the lead-lag filter
            T13_ll = 0.0 # time parameter for the lead-lag filter
            T14_ll = 0.0 # time parameter for the lead-lag filter
            T15_ll = 0.0 # time parameter for the lead-lag filter
            T16_ll = 0.0 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt = lead_lag_8th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                x8,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
                T15_ll,
                T16_ll,
            )
        end
    elseif N == 5
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 5.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 10.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 10.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 5.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = T1^5 # time parameter for the lead-lag filter
            T6_ll = 5.0 * T2 # time parameter for the lead-lag filter
            T7_ll = 10.0 * T2^2 # time parameter for the lead-lag filter
            T8_ll = 10.0 * T2^3 # time parameter for the lead-lag filter
            T9_ll = 5.0 * T2^4 # time parameter for the lead-lag filter
            T10_ll = T2^5 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt = lead_lag_5th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
            )
        end
    elseif N == 6
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 6.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 15.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 20.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 15.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 6.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = T1^6 # time parameter for the lead-lag filter
            T7_ll = 6.0 * T2 # time parameter for the lead-lag filter
            T8_ll = 15.0 * T2^2 # time parameter for the lead-lag filter
            T9_ll = 20.0 * T2^3 # time parameter for the lead-lag filter
            T10_ll = 15.0 * T2^4 # time parameter for the lead-lag filter
            T11_ll = 6.0 * T2^5 # time parameter for the lead-lag filter
            T12_ll = T2^6 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt = lead_lag_6th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
            )
        end
    elseif N == 7
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 7.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 21.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 35.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 35.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 21.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 7.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = T1^7 # time parameter for the lead-lag filter
            T8_ll = 7.0 * T2 # time parameter for the lead-lag filter
            T9_ll = 21.0 * T2^2 # time parameter for the lead-lag filter
            T10_ll = 35.0 * T2^3 # time parameter for the lead-lag filter
            T11_ll = 35.0 * T2^4 # time parameter for the lead-lag filter
            T12_ll = 21.0 * T2^5 # time parameter for the lead-lag filter
            T13_ll = 7.0 * T2^6 # time parameter for the lead-lag filter
            T14_ll = T2^7 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt = lead_lag_7th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
            )
        end
    elseif N == 8
        if M == 0
            y = u
        elseif M == 1
            T1_ll = 8.0 * T1 # time parameter for the lead-lag filter
            T2_ll = 28.0 * T1^2 # time parameter for the lead-lag filter
            T3_ll = 56.0 * T1^3 # time parameter for the lead-lag filter
            T4_ll = 70.0 * T1^4 # time parameter for the lead-lag filter
            T5_ll = 56.0 * T1^5 # time parameter for the lead-lag filter
            T6_ll = 28.0 * T1^6 # time parameter for the lead-lag filter
            T7_ll = 8.0 * T1^7 # time parameter for the lead-lag filter
            T8_ll = T1^8 # time parameter for the lead-lag filter
            T9_ll = 8.0 * T2 # time parameter for the lead-lag filter
            T10_ll = 28.0 * T2^2 # time parameter for the lead-lag filter
            T11_ll = 56.0 * T2^3 # time parameter for the lead-lag filter
            T12_ll = 70.0 * T2^4 # time parameter for the lead-lag filter
            T13_ll = 56.0 * T2^5 # time parameter for the lead-lag filter
            T14_ll = 28.0 * T2^6 # time parameter for the lead-lag filter
            T15_ll = 8.0 * T2^7 # time parameter for the lead-lag filter
            T16_ll = T2^8 # time parameter for the lead-lag filter

            y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt = lead_lag_8th(
                u,
                x1,
                x2,
                x3,
                x4,
                x5,
                x6,
                x7,
                x8,
                T1_ll,
                T2_ll,
                T3_ll,
                T4_ll,
                T5_ll,
                T6_ll,
                T7_ll,
                T8_ll,
                T9_ll,
                T10_ll,
                T11_ll,
                T12_ll,
                T13_ll,
                T14_ll,
                T15_ll,
                T16_ll,
            )
        end
    end

    return y, dx1dt, dx2dt, dx3dt, dx4dt, dx5dt, dx6dt, dx7dt, dx8dt
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
    y::Z,
    K::ACCEPTED_REAL_TYPES,
    ::Float64,
    y_min::ACCEPTED_REAL_TYPES,
    y_max::ACCEPTED_REAL_TYPES,
) where {Z <: ACCEPTED_REAL_TYPES}
    dydt_scaled = K * u
    y_sat = clamp(y, y_min, y_max)
    return y_sat, dydt_scaled
end

# Does not accept T = 0
function integrator_windup(
    u::Z,
    y::Z,
    K::ACCEPTED_REAL_TYPES,
    T::ACCEPTED_REAL_TYPES,
    y_min::ACCEPTED_REAL_TYPES,
    y_max::ACCEPTED_REAL_TYPES,
) where {Z <: ACCEPTED_REAL_TYPES}
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
    y::Z,
    K::Float64,
    ::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    dydt_scaled = K * u
    y_sat = clamp(y, y_min, y_max)
    upper_lim = (y >= y_max) && (dydt_scaled > 0)
    lower_lim = (y <= y_min) && (dydt_scaled < 0)
    binary_logic = upper_lim || lower_lim ? 0.0 : 1.0
    return y_sat, binary_logic * dydt_scaled
end

# Does not accept T = 0
function integrator_nonwindup(
    u::Z,
    y::Z,
    K::Float64,
    T::Float64,
    y_min::Float64,
    y_max::Float64,
) where {Z <: ACCEPTED_REAL_TYPES}
    y_sat, dydt_scaled = integrator_nonwindup_mass_matrix(u, y, K, T, y_min, y_max)
    return y_sat, (1.0 / T) * dydt_scaled
end
