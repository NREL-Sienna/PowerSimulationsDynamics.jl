abstract type TurbineGov <: GeneratorComponent end

@def TG_ports begin
      state_input = [:ω]
      inner_input = Vector{Int64}()
end



"""
Parameters of a fixed Turbine Governor that returns a fixed mechanical torque
given by the product of P_ref*efficiency
"""


mutable struct TGFixed <: TurbineGov
      efficiency::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function TGFixed(efficiency::Float64)
            n_states = 0
            states = Vector{Symbol}()
            @TG_ports
            new(efficiency, n_states, states, Ports(state_input, inner_input))
      end
end

"""
Parameters of a Turbine Governor Type I.

# Conmutable structor
```julia
TGTypeI(R, Ts, Tc, T3, T4, T5, P_min, P_max)
```

#Arguments
* `R`::Float64 : Droop parameter
* `Ts`::Float64 : Governor time constant
* `Tc`::Float64 : Servo time constant
* `T3`::Float64 : Transient gain time constant
* `T4`::Float64 : Power fraction time constant
* `T5`::Float64 : Reheat time constant
* `P_min`::Float64 : Min Power into the Governor
* `P_max`::Float64 : Max Power into the Governor

"""


mutable struct TGTypeI <: TurbineGov
      R::Float64
      Ts::Float64
      Tc::Float64
      T3::Float64
      T4::Float64
      T5::Float64
      P_min::Float64
      P_max::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function TGTypeI(R::Float64,
                       Ts::Float64,
                       Tc::Float64,
                       T3::Float64,
                       T4::Float64,
                       T5::Float64,
                       P_min::Float64,
                       P_max::Float64)

            n_states = 3
            states = [:x_g1, :x_g2, :x_g3]
            @TG_ports
            new(R, Ts, Tc, T3, T4, T5, P_min,
                P_max, n_states, states, Ports(state_input, inner_input))
      end
end



"""
Parameters of a Turbine Governor Type II.

# Conmutable structor
```julia
TGTypeI(R, T1, T2, τ_min, τ_max)
```

#Arguments
* `R`::Float64 : Droop parameter
* `T1`::Float64 : Transient gain time constant
* `T2`::Float64 : Governor time constant
* `τ_min`::Float64 : Min turbine torque output
* `τ_max`::Float64 : Max turbine torque output

"""


mutable struct TGTypeII <: TurbineGov
      R::Float64
      T1::Float64
      T2::Float64
      τ_min::Float64
      τ_max::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function TGTypeII(R::Float64,
                       T1::Float64,
                       T2::Float64,
                       τ_min::Float64,
                       τ_max::Float64)

            n_states = 1
            states = [:xg]
            @TG_ports
            new(R, T1, T2, τ_min, τ_max, n_states, states, Ports(state_input, inner_input))
      end
end
