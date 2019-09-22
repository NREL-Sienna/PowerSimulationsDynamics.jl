abstract type AVR <: GeneratorComponent end

@def avr_ports begin
      state_input = Vector{Symbol}()
      inner_input = [V_pss_var, VI_gen_var, VR_gen_var]
end

"""
Parameters of a AVR that returns a fixed voltage
"""

mutable struct AVRFixed <: AVR
      Emf::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function AVRFixed(Emf::Float64)
            n_states = 0
            states = Vector{Symbol}()
            @avr_ports
            new(Emf, n_states, states, Ports(state_input, inner_input))
      end
end


"""
Parameters of a simple proportional AVR in the derivative of EMF
i.e. an integrator controller on EMF
"""
mutable struct AVRSimple <: AVR
     Kv::Float64
     n_states::Int64
     states::Vector{Symbol}
     ports::Ports

     function AVRSimple(Kv::Float64)
        n_states = 1
        states = [:Vf]
        @avr_ports
        new(Kv, n_states, states, Ports(state_input, inner_input))
     end

end

"""
Parameters of an Automatic Voltage Regulator Type I - Resembles IEEE Type DC1

# Conmutable structor
```julia
AVRTypeI(Ka, Ke, Kf, Ta, Tf, Te, Tr, Vr_max, Vr_min, Ae, Be)
```

#Arguments
* `Ka`::Float64 : Amplifier Gain
* `Ke`::Float64 : Field circuit integral deviation
* `Kf`::Float64 : Stabilizer Gain in s * pu/pu
* `Ta`::Float64 : Amplifier Time Constant in s
* `Te`::Float64 : Field Circuit Time Constant in s
* `Tf`::Float64 : Stabilizer Time Constant in s
* `Tr`::Float64 : Voltage Measurement Time Constant in s
* `Vr_max`::Float64 : Maximum regulator voltage in pu
* `Vr_min`::Float64 : Minimum regulator voltage in pu
* `Ae`::Float64 : 1st ceiling coefficient
* `Be`::Float64 : 2nd ceiling coefficient
"""


mutable struct AVRTypeI <: AVR
      Ka::Float64
      Ke::Float64
      Kf::Float64
      Ta::Float64
      Te::Float64
      Tf::Float64
      Tr::Float64
      Vr_max::Float64
      Vr_min::Float64
      Ae::Float64
      Be::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function AVRTypeI(Ka::Float64,
                        Ke::Float64,
                        Kf::Float64,
                        Ta::Float64,
                        Te::Float64,
                        Tf::Float64,
                        Tr::Float64,
                        Vr_max::Float64,
                        Vr_min::Float64,
                        Ae::Float64,
                        Be::Float64) where {P<:PSS}

         n_states = 4
         states = [:Vf, :Vr1, :Vr2, :Vm]
         @avr_ports
         new(Ka, Ke, Kf, Ta, Te, Tf, Tr, Vr_max,
             Vr_min, Ae, Be, n_states, states, Ports(state_input, inner_input))
      end
end


"""
Parameters of an Automatic Voltage Regulator Type II -  Typical static exciter model

# Conmutable structor
```julia
AVRTypeII(K0, T1, T2, T3, T4, Te, Tr, Vr_max, Vr_min, Ae, Be)
```

#Arguments
* `K0`::Float64 : Regulator Gain
* `T1`::Float64 : First Pole in s
* `T2`::Float64 : First zero in s
* `T3`::Float64 : First Pole in s
* `T4`::Float64 : First zero in s
* `Te`::Float64 : Field Circuit Time Constant in s
* `Tr`::Float64 : Voltage Measurement Time Constant in s
* `Vr_max`::Float64 : Maximum regulator voltage in pu
* `Vr_min`::Float64 : Minimum regulator voltage in pu
* `Ae`::Float64 : 1st ceiling coefficient
* `Be`::Float64 : 2nd ceiling coefficient
"""

mutable struct AVRTypeII <: AVR
      K0::Float64
      T1::Float64
      T2::Float64
      T3::Float64
      T4::Float64
      Te::Float64
      Tr::Float64
      Vr_max::Float64
      Vr_min::Float64
      Ae::Float64
      Be::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function AVRTypeII(K0::Float64,
                        T1::Float64,
                        T2::Float64,
                        T3::Float64,
                        T4::Float64,
                        Te::Float64,
                        Tr::Float64,
                        Vr_max::Float64,
                        Vr_min::Float64,
                        Ae::Float64,
                        Be::Float64) where {P<:PSS}

         n_states = 4
         states = [:Vf, :Vr1, :Vr2, :Vm]
         @avr_ports
         new(K0, T1, T2, T3, T4, Te, Tr, Vr_max,
             Vr_min, Ae, Be, n_states, states, Ports(state_input, inner_input))
      end
end



"""
Parameters of an Automatic Voltage Regulator Type II from PSAT Manual -  Typical static exciter model

# Conmutable structor
```julia
AVRTypeII(K0, T1, T2, T3, T4, Te, Tr, Vr_max, Vr_min, Ae, Be)
```

#Arguments
* `Kf`::Float64 : Regulator Gain
* `Ta`::Float64 : Amplifier time constant
* `Kf`::Float64 : Stabilizer gain
* `Tf`::Float64 : Stabilizer time constant
* `Te`::Float64 : Field Circuit Time Constant in s
* `Tr`::Float64 : Voltage Measurement Time Constant in s
* `Vr_max`::Float64 : Maximum regulator voltage in pu
* `Vr_min`::Float64 : Minimum regulator voltage in pu
* `Ae`::Float64 : 1st ceiling coefficient
* `Be`::Float64 : 2nd ceiling coefficient
"""

mutable struct AVRTypeIIManual <: AVR
      Ka::Float64
      Ta::Float64
      Kf::Float64
      Tf::Float64
      Te::Float64
      Tr::Float64
      Vr_max::Float64
      Vr_min::Float64
      Ae::Float64
      Be::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function AVRTypeIIManual(Ka::Float64,
                        Ta::Float64,
                        Kf::Float64,
                        Tf::Float64,
                        Te::Float64,
                        Tr::Float64,
                        Vr_max::Float64,
                        Vr_min::Float64,
                        Ae::Float64,
                        Be::Float64) where {P<:PSS}

         n_states = 4
         states = [:Vf, :Vr1, :Vr2, :Vm]
         @avr_ports
         new(Ka, Ta, Kf, Tf, Te, Tr, Vr_max,
             Vr_min, Ae, Be, n_states, states, Ports(state_input, inner_input))
      end
end
