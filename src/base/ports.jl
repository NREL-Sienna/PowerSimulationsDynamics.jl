function Ports(states::Vector{Symbol}, inner::Vector)
    return Dict(:states=>states, :inner=>Int.(inner))
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
    state_input = Vector{Symbol}()
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
    inner_input = [md_var, mq_var, Vdc_var, Vr_cnv_var, Vi_cnv_var]
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
        Vr_cnv_var,
        Vi_cnv_var,
        θ_oc_var,
        Vr_filter_var,
        Vi_filter_var,
    ]
    return Ports(state_input, inner_input)
end

#### Freq. Estimator Ports ####

function Ports(::PSY.FrequencyEstimator)
    state_input = [:vr_filter, :vi_filter, :θ_oc]
    #TODO: Move PLL to PCC, i.e. move v_cap (D'Arco v_o), to inner inputs
    inner_input = [Vr_filter_var, Vi_filter_var, θ_oc_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end

#### Outer Control Ports ####
function Ports(::PSY.OuterControl)
    state_input = [:vd_pll, :vq_pll, :ε_pll, :vr_filter, :vi_filter, :ir_filter, :ii_filter]
    inner_input = [Vr_filter_var, Vi_filter_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end

#### Inner Control Ports ####
function Ports(::PSY.InnerControl)
    state_input = [:ir_filter, :ii_filter, :ir_cnv, :ii_cnv, :vr_filter, :vi_filter]
    inner_input = [Vr_filter_var, Vi_filter_var, V_oc_var, ω_oc_var]
    return Ports(state_input, inner_input)
end
