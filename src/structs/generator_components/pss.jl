abstract type PSS <: GeneratorComponent end

@def pss_ports begin
      state_input = [:ω]
      inner_input = [τe_var, VR_gen_var]
  end


"""
Parameters of a PSS that returns a fixed voltage to add to
the reference for the AVR

# Conmutable structor
```julia
PSSFixed(Vs)
```

# Arguments
* `Vs`::Float64 : Fixed voltage stabilization signal

"""
mutable struct PSSFixed <: PSS
      Vs::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function PSSFixed(Vs::Float64)
            n_states = 0
            states = Vector{Symbol}()
            @pss_ports
            new(Vs, n_states, states, Ports(state_input, inner_input))
      end
end


"""
Parameters of a PSS that returns a proportional droop voltage to add to
the reference for the AVR

# Conmutable structor
```julia
PSSSimple(K_ω, K_p)
```

# Arguments
* `K_ω`::Float64 : Proportional gain for frequency
* `K_p`::Float64 : Proportional gain for active power

"""
mutable struct PSSSimple <: PSS
      K_ω::Float64
      K_p::Float64
      n_states::Int64
      states::Vector{Symbol}
      ports::Ports

      function PSSSimple(K_ω::Float64, K_p::Float64)
            n_states = 0
            states = Vector{Symbol}()
            @pss_ports
            new(K_ω, K_p, n_states, states, Ports(state_input, inner_input))
      end
end
