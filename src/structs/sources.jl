"""
This struct acts as an infinity bus.
"""
struct StaticSource <: PSY.Injection
    number::Int64
    name::Symbol
    bus::PSY.Bus
    V_R::Float64
    V_I::Float64
    X_th::Float64
end
