### Port struct ###
mutable struct Ports
    state::Vector{Symbol}
    inner::Vector
end

#### AVR Ports ####
function Ports(::PSY.AVR)
    state_input = Vector{Symbol}()
    inner_input = [V_pss_var, VI_gen_var, VR_gen_var]
    return Ports(state_input, inner_input)
end

#### Machine Ports ####
function Ports(::PSY.Machine)
    state_input = [:δ, :ω, :Vf]
    inner_input = [VI_gen_var, VR_gen_var]
    return Ports(state_input, inner_input)
end

#### PSS Ports ####
function Ports(::PSY.PSS)
    state_input = [:ω]
    inner_input = [τe_var, VR_gen_var]
    return Ports(state_input, inner_input)
end

#### Shaft Ports ####
function Ports(::PSY.Shaft)
    state_input = Vector{Int64}()
    inner_input = [τe_var, τm_var]
    return Ports(state_input, inner_input)
end

#### Turbine Governor Ports ####
function Ports(::PSY.TurbineGov)
    state_input = [:ω]
    inner_input = Vector{Int64}()
    return Ports(state_input, inner_input)
end

##################################################
############ Inverter Components ################

#### Converter Ports ####
function Ports(::PSY.Converter)
    state_input = Vector{Symbol}()
    inner_input = [md_var, mq_var, Vdc_var, Vd_cnv_var, Vq_cnv_var]
    return Ports(state_input, inner_input)
end

#### DC Source Ports ####
function Ports(::PSY.DCSource)
    state_input = Vector{Symbol}()
    inner_input = Vector{Int64}()
    return Ports(state_input, inner_input)
end

#### Filter Ports ####
function Ports(::PSY.Filter)
    #TODO: If converter has dynamics, need to connect state_input
    state_input = [:θ_oc] #[:Vd_c, :Vq_c] #, :Id_c, :Iq_c]
    inner_input = [
        VR_inv_var,
        VI_inv_var,
        Vd_cnv_var,
        Vq_cnv_var,
        θ_oc_var,
        Vd_filter_var,
        Vq_filter_var,
    ]
    return Ports(state_input, inner_input)
end

#### Freq. Estimator Ports ####

function Ports(::PSY.FrequencyEstimator)
    state_input = [:vd_filter, :vq_filter, :θ_oc]
    #TODO: Move PLL to PCC, i.e. move v_cap (D'Arco v_o), to inner inputs
    inner_input = [Vd_filter_var, Vq_filter_var, θ_oc_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end

#### Outer Control Ports ####
function Ports(::PSY.OuterControl)
    state_input = [:vd_pll, :vq_pll, :ε_pll, :vd_filter, :vq_filter, :id_filter, :iq_filter]
    inner_input = [Vd_filter_var, Vq_filter_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end

#### Inner Control Ports ####
function Ports(::PSY.InnerControl)
    state_input = [:id_filter, :iq_filter, :id_cnv, :iq_cnv, :vd_filter, :vq_filter]
    inner_input = [Vd_filter_var, Vq_filter_var, V_oc_var, ω_oc_var]
    return Ports(state_input, inner_input)
end
