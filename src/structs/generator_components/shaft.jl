abstract type Shaft <: GeneratorComponent end

@def shaft_ports begin
    state_input = Vector{Int64}()
    inner_input = [τe_var, τm_var]
end


"""
Parameters of single mass shaft model. Typically represents the rotor mass.

# Conmutable structor
```julia
SingleMass(H, D)
```

# Arguments
* `H`   ::Float64 : Rotor inertia constant in MWs/MVA
* `D`   ::Float64 : Rotor natural damping in pu

"""
mutable struct SingleMass <: Shaft
    H::Float64
    D::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

    function SingleMass(H::Float64,
                        D::Float64)

        n_states = 2
        states = [:δ, :ω]
        @shaft_ports

        new(H,
            D,
            n_states,
            states,
            Ports(state_input, inner_input))
    end

end

"""
Parameters of 5 mass-spring shaft model.
It contains a High-Pressure (HP) steam turbine, Intermediate-Pressure (IP)
steam turbine, Low-Pressure (LP) steam turbine, the Rotor and
an Exciter (EX) mover.

# Conmutable structor
```julia
FiveMassShaft(H, H_hp, H_ip, H_lp, H_ex,
              D, D_hp, D_ip, D_lp, D_ex,
              D_12, D_23, D_34, D_45,
              K_hp, K_ip, K_lp, K_ex)
```

# Arguments
* `H`   ::Float64 : Rotor inertia constant in MWs/MVA
* `H_hp`::Float64 : High pressure turbine inertia constant in MWs/MVA
* `H_ip`::Float64 : Intermediate pressure turbine inertia constant in MWs/MVA
* `H_lp`::Float64 : Low pressure turbine inertia constant in MWs/MVA
* `H_ex`::Float64 : Exciter inertia constant in MWs/MVA
* `D`   ::Float64 : Rotor natural damping in pu
* `D_hp`::Float64 : High pressure turbine natural damping in pu
* `D_ip`::Float64 : Intermediate pressure turbine natural damping in pu
* `D_lp`::Float64 : Low pressure turbine natural damping in pu
* `D_ex`::Float64 : Exciter natural damping in pu
* `D_12`::Float64 : High-Intermediate pressure turbine damping
* `D_23`::Float64 : Intermediate-Low pressure turbine damping
* `D_34`::Float64 : Low pressure turbine-Rotor damping
* `D_45`::Float64 : Rotor-Exciter damping
* `K_hp`::Float64 : High pressure turbine angle coefficient
* `K_ip`::Float64 : Intermediate pressure turbine angle coefficient
* `K_lp`::Float64 : Low pressure turbine angle coefficient
* `K_ex`::Float64 : Exciter angle coefficient

"""
mutable struct FiveMassShaft <: Shaft
    H::Float64
    H_hp::Float64
    H_ip::Float64
    H_lp::Float64
    H_ex::Float64
    D::Float64
    D_hp::Float64
    D_ip::Float64
    D_lp::Float64
    D_ex::Float64
    D_12::Float64
    D_23::Float64
    D_34::Float64
    D_45::Float64
    K_hp::Float64
    K_ip::Float64
    K_lp::Float64
    K_ex::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

    function FiveMassShaft(H::Float64,
                           H_hp::Float64,
                           H_ip::Float64,
                           H_lp::Float64,
                           H_ex::Float64,
                           D::Float64,
                           D_hp::Float64,
                           D_ip::Float64,
                           D_lp::Float64,
                           D_ex::Float64,
                           D_12::Float64,
                           D_23::Float64,
                           D_34::Float64,
                           D_45::Float64,
                           K_hp::Float64,
                           K_ip::Float64,
                           K_lp::Float64,
                           K_ex::Float64)

        n_states = 10
        states = [:δ, :ω,
                  :δ_hp, :ω_hp,
                  :δ_ip, :ω_ip,
                  :δ_lp, :ω_lp,
                  :δ_ex, :ω_ex]
        @shaft_ports

        new(H,
            H_hp,
            H_ip,
            H_lp,
            H_ex,
            D,
            D_hp,
            D_ip,
            D_lp,
            D_ex,
            D_12,
            D_23,
            D_34,
            D_45,
            K_hp,
            K_ip,
            K_lp,
            K_ex,
            n_states,
            states,
            Ports(state_input, inner_input))
    end

end
