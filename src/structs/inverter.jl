@enum inverter_inner_vars begin
    md_var = 1
    mq_var = 2
    Vdc_var = 3
    Vdo_var = 5
    Vqo_var = 6
    ω_freq_estimator_var = 7
    v_control_var = 8
    ω_control_var = 9
    δdqRI_var = 10 # alias for δθ_vsm, TODO: should δθ_vsm ever be accessed as external state? RH: Most likely
    VR_inv_var = 11
    VI_inv_var = 12
    Vdcnv_var = 13
    Vqcnv_var = 14
end

Base.to_index(ix::inverter_inner_vars) = Int(ix)

include("inverter_components/outer_control.jl")
include("inverter_components/voltage_source_control.jl")
include("inverter_components/converter.jl")
include("inverter_components/dc_source.jl")
include("inverter_components/frequency_estimator.jl")
include("inverter_components/filter.jl")

mutable struct DynInverter{C <: Converter,
                           O <: OuterControl,
                           VC <: VSControl,
                           DC<: DCSource,
                           P <: FrequencyEstimator,
                           F <: Filter} <: DynInjection
    number::Int64
    name::Symbol
    bus::PSY.Bus
    ω_ref::Float64
    V_ref::Float64
    P_ref::Float64
    Q_ref::Float64
    MVABase::Float64
    converter::C
    outercontrol::O
    vscontrol::VC
    dc_source::DC
    freq_estimator::P
    filter::F #add MVAbase here
    n_states::Int64
    states::Vector{Symbol}
    inner_vars::Vector{Float64}
    local_state_ix::Dict{InverterComponent,Vector{Int64}}
    input_port_mapping::Dict{InverterComponent,Vector{Int64}}

        function DynInverter(number::Int64,
                            name::Symbol,
                            bus::PSY.Bus,
                            ω_ref::Float64,
                            V_ref::Float64,
                            P_ref::Float64,
                            Q_ref::Float64,
                            MVABase::Float64,
                            converter::C,
                            outercontrol::O,
                            vscontrol::VC,
                            dc_source::DC,
                            freq_estimator::P,
                            filter::F) where {C<:Converter,
                                              O<:OuterControl,
                                              VC<:VSControl,
                                              DC<:DCSource,
                                              P<:FrequencyEstimator,
                                              F<:Filter}

            n_states = (converter.n_states +
                       outercontrol.n_states +
                       vscontrol.n_states +
                       dc_source.n_states +
                       freq_estimator.n_states +
                       filter.n_states)

            states = vcat(converter.states,
                        outercontrol.states,
                        vscontrol.states,
                        dc_source.states,
                        freq_estimator.states,
                        filter.states)

            components = [converter,
                          outercontrol,
                          vscontrol,
                          dc_source,
                          freq_estimator,
                          filter]

            new{C, O, VC, DC, P, F}(number,
                                    name,
                                    bus,
                                    ω_ref,
                                    V_ref,
                                    P_ref,
                                    Q_ref,
                                    MVABase,
                                    converter,
                                    outercontrol,
                                    vscontrol,
                                    dc_source,
                                    freq_estimator,
                                    filter,
                                    n_states,
                                    states,
                                    zeros(14),
                                    _make_state_mapping(components, states),
                                    _make_port_mapping(components, states))

        end
end

get_inverter_Sbase(device::DynInverter) = device.converter.s_rated
get_inverter_Vref(device::DynInverter) = device.V_ref
