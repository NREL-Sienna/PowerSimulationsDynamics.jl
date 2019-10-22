abstract type OuterControl <: InverterComponent end

abstract type ActivePowerControl end
abstract type ReativePowerControl end

@def outercontrol_ports begin
    state_input = [:vpll_d, :vpll_q, :ε_pll, :vd_cap, :vq_cap, :id_o, :iq_o]
    inner_input = [Vdo_var, Vdo_var, ω_freq_estimator_var]
end

"""
Parameters of a Outer-Loop controller using a virtual inertia with VSM for active power controller
and a reactive power droop controller.
# Conmutable structor
```julia
VirtualInertiaQDroop(A, R)
```

# Arguments
*  `A`::Float64 : Active power controller using virtual inertia with VSM
*  `R`::Float64 : Reactive power controller using reactive power droop
"""
mutable struct VirtualInertiaQdroop{A <: ActivePowerControl,
                                   R <: ReativePowerControl} <: OuterControl
    active_power::A
    reactive_power::R
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function VirtualInertiaQdroop(active_power::A,
                                     reactive_power::R) where {A <: ActivePowerControl,
                                                               R <: ReativePowerControl}

            total_states = active_power.n_states + reactive_power.n_states
            @outercontrol_ports

            new{A,R}(active_power,
                     reactive_power,
                     total_states,
                     vcat(active_power.states, reactive_power.states),
                     Ports(state_input, inner_input))

        end

end

"""
Parameters of a Virtual Inertia with SRF using VSM for active power controller
# Conmutable structor
```julia
VirtualInertia(Ta, kd, kw, ωb)
```

# Arguments
*  `Ta`::Float64 : VSM inertia constant
*  `kd`::Float64 : VSM damping constant
*  `kw`::Float64 : frequency droop gain
*  `ωb`::Float64 : rated angular frequency
"""
mutable struct VirtualInertia <: ActivePowerControl
    Ta::Float64
    kd::Float64
    kω::Float64
    ωb::Float64
    n_states::Int64
    states::Vector{Symbol}

        function VirtualInertia(Ta::Float64,
                                kd::Float64,
                                kω::Float64,
                                ωb::Float64)

            new(Ta,
                kd,
                kω,
                ωb,
                2,
                [:δω_vsm, :δθ_vsm])

        end

end

"""
Parameters of a Reactive Power droop controller
# Conmutable structor
```julia
ReactivePowerDroop(kq, ωf)
```

# Arguments
*  `kq`::Float64 : reactive power droop gain
*  `ωf`::Float64 : reactive power filter cutoff frequency (rad/sec)
"""
mutable struct ReactivePowerDroop <: ReativePowerControl
    kq::Float64
    ωf::Float64
    n_states::Int64
    states::Vector{Symbol}

        function ReactivePowerDroop(kq::Float64,
                                   ωf::Float64)

            new(kq, ωf, 1, [:qm])
        end
end
