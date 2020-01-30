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
    inner_input = [md_var, mq_var, Vdc_var, Vdcnv_var, Vqcnv_var]
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
    state_input = [:δθ_vsm] #[:Vd_c, :Vq_c] #, :Id_c, :Iq_c]
    inner_input =
        [VR_inv_var, VI_inv_var, Vdcnv_var, Vqcnv_var, δdqRI_var, Vdo_var, Vqo_var]
    return Ports(state_input, inner_input)
end


#### Freq. Estimator Ports ####

function Ports(::PSY.FrequencyEstimator)
    state_input = [:vd_cap, :vq_cap, :δθ_vsm]
    #TODO: Move PLL to PCC, i.e. move v_cap (D'Arco v_o), to inner inputs
    inner_input = [Vdo_var, Vqo_var, δdqRI_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end

#### Outer Control Ports ####
function Ports(::PSY.OuterControl)
    state_input = [:vpll_d, :vpll_q, :ε_pll, :vd_cap, :vq_cap, :id_o, :iq_o]
    inner_input = [Vdo_var, Vdo_var, ω_freq_estimator_var]
    return Ports(state_input, inner_input)
end


#### Inner Control Ports ####
function Ports(::PSY.VSControl)
    state_input = [:id_o, :iq_o, :id_c, :iq_c, :vd_cap, :vq_cap]
    inner_input = [Vdo_var, Vqo_var, v_control_var, ω_control_var]
    return Ports(state_input, inner_input)
end
