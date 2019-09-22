abstract type FrequencyEstimator <: InverterComponent end

@def freq_estimation_ports begin
    state_input = [:vd_cap,:vq_cap, :δθ_vsm]
    #TODO: Move PLL to PCC, i.e. move v_cap (D'Arco v_o), to inner inputs
    inner_input = [Vdo_var,Vqo_var,δdqRI_var,ω_freq_estimator_var]
end

"""
*  `ω_lp`::Float64 : PLL low-pass filter frequency (rad/sec)
*  `kp_pll`::Float64 : PLL proportional gain
*  `ki_pll`::Float64 : PLL integral gain
"""
mutable struct PLL <: FrequencyEstimator
    ω_lp::Float64
    kp_pll::Float64
    ki_pll::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function PLL(ω_lp::Float64,
                     kp_pll::Float64,
                     ki_pll::Float64)

            @freq_estimation_ports

            new(ω_lp,
                kp_pll,
                ki_pll,
                4,
                [:vpll_d, :vpll_q, :ε_pll, :δθ_pll],
                Ports(state_input, inner_input))

        end

end
