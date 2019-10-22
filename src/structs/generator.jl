@enum generator_inner_vars begin
    τe_var = 1
    τm_var = 2
    Vf_var = 3
    V_pss_var = 4
    VR_gen_var = 5
    VI_gen_var = 6
    ψd_var = 7
    ψq_var = 8
end

Base.to_index(ix::generator_inner_vars) = Int(ix)
include("generator_components/machine.jl")
include("generator_components/shaft.jl")
include("generator_components/pss.jl")
include("generator_components/avr.jl")
include("generator_components/turbinegov.jl")

"""
"""
mutable struct DynGenerator{M <: Machine,
                            S <: Shaft,
                            A <: AVR,
                            TG<: TurbineGov,
                            P <: PSS} <: DynInjection
    number::Int64
    name::Symbol
    bus::PSY.Bus
    ω_ref::Float64
    V_ref::Float64
    P_ref::Float64
    machine::M
    shaft::S
    avr::A
    tg::TG
    pss::P
    n_states::Int64
    states::Vector{Symbol}
    inner_vars::Vector{Float64}
    local_state_ix::Dict{GeneratorComponent,Vector{Int64}}
    input_port_mapping::Dict{GeneratorComponent,Vector{Int64}}

        function DynGenerator(number::Int64,
                              name::Symbol,
                              bus::PSY.Bus,
                              ω_ref::Float64,
                              V_ref::Float64,
                              P_ref::Float64,
                              machine::M,
                              shaft::S,
                              avr::A,
                              tg::TG,
                              pss::P) where {M <: Machine,
                                             S <: Shaft,
                                             A <: AVR,
                                             TG <: TurbineGov,
                                             P <: PSS}

            n_states = (machine.n_states +
                       shaft.n_states +
                       avr.n_states +
                       tg.n_states +
                       pss.n_states)

            states = vcat(machine.states,
                          shaft.states,
                          avr.states,
                          tg.states,
                          pss.states)

            components = [machine, shaft, avr, tg, pss]

            new{M,S,A,TG,P}(number,
                        name,
                        bus,
                        ω_ref,
                        V_ref,
                        P_ref,
                        machine,
                        shaft,
                        avr,
                        tg,
                        pss,
                        n_states,
                        states,
                        zeros(8),
                        _make_state_mapping(components, states),
                        _make_port_mapping(components, states))

        end
  end
