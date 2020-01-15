### Port struct ###
mutable struct Ports <: DynDevice
    state::Vector{Symbol}
    inner::Vector
end

##################################################
############### Port macros ######################
##################################################


##################################################
############ Generator Components ################

#### AVR Ports ####

@def avr_ports begin
    state_input = Vector{Symbol}()
    inner_input = [V_pss_var, VI_gen_var, VR_gen_var]
end

#### Machine Ports ####

@def machine_ports begin
    state_input = [:δ, :ω, :Vf]
    inner_input = [VI_gen_var, VR_gen_var]
end

#### PSS Ports ####

@def pss_ports begin
    state_input = [:ω]
    inner_input = [τe_var, VR_gen_var]
end

#### Shaft Ports ####

@def shaft_ports begin
    state_input = Vector{Int64}()
    inner_input = [τe_var, τm_var]
end

#### Turbine Governor Ports ####

@def TG_ports begin
    state_input = [:ω]
    inner_input = Vector{Int64}()
end


##################################################
############ Inverter Components ################

#### Converter Ports ####

@def converter_ports begin
    state_input = Vector{Symbol}()
    inner_input = [md_var, mq_var, Vdc_var, Vdcnv_var,Vqcnv_var]
end

#### DC Source Ports ####

@def dcsource_ports begin
    state_input = Vector{Symbol}()
    inner_input = Vector{Int64}()
end


#### Filter Ports ####

@def filter_ports begin
    #TODO: If converter has dynamics, need to connect state_input
    state_input = [:δθ_vsm] #[:Vd_c, :Vq_c] #, :Id_c, :Iq_c]
    inner_input = [VR_inv_var,VI_inv_var,Vdcnv_var,Vqcnv_var,δdqRI_var,Vdo_var,Vqo_var]
end


#### Freq. Estimator Ports ####

@def freq_estimation_ports begin
    state_input = [:vd_cap,:vq_cap, :δθ_vsm]
    #TODO: Move PLL to PCC, i.e. move v_cap (D'Arco v_o), to inner inputs
    inner_input = [Vdo_var,Vqo_var,δdqRI_var,ω_freq_estimator_var]
end

#### Outer Control Ports ####

@def outercontrol_ports begin
    state_input = [:vpll_d, :vpll_q, :ε_pll, :vd_cap, :vq_cap, :id_o, :iq_o]
    inner_input = [Vdo_var, Vdo_var, ω_freq_estimator_var]
end


#### Inner Control Ports ####

@def vscontrol_ports begin
    state_input = [:id_o, :iq_o, :id_c, :iq_c, :vd_cap, :vq_cap]
    inner_input = [Vdo_var, Vqo_var, v_control_var, ω_control_var]
end
