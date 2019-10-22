abstract type VSControl <: InverterComponent end

@def vscontrol_ports begin
    state_input = [:id_o, :iq_o, :id_c, :iq_c, :vd_cap, :vq_cap]
    inner_input = [Vdo_var, Vqo_var, v_control_var, ω_control_var]
end

"""
Parameters of an inner loop controller using virtual impedance, voltage controller and current controller.
# Conmutable structor
```julia
CombinedVIwithVZ(kpv, kiv, kffv, rv, lv, kpc, kic, kffi, ωad, kad)
```

# Arguments
*  `kpv`::Float64 : voltage controller proportional gain
*  `kiv`::Float64 : voltage controller integral gain
*  `kffv`::Float64 : Binary variable to enable feed-forward gain of voltage.
*  `rv`::Float64 : virtual resistance
*  `lv`::Float64 : virtual inductance
*  `kpc`::Float64 : current controller proportional gain
*  `kic`::Float64 : current controller integral gain
*  `kffi`::Float64 : Binary variable to enable feed-forward gain of current
*  `ωad`::Float64 : active damping filter cutoff frequency (rad/sec)
*  `kad`::Float64 : active damping gain
"""
mutable struct CombinedVIwithVZ <: VSControl #(D'Arco EPSR122 Model)
    kpv::Float64
    kiv::Float64
    kffv::Float64
    rv::Float64
    lv::Float64
    kpc::Float64
    kic::Float64
    kffi::Float64
    ωad::Float64
    kad::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function CombinedVIwithVZ(kpv::Float64,
                                  kiv::Float64,
                                  kffv::Float64,
                                  rv::Float64,
                                  lv::Float64,
                                  kpc::Float64,
                                  kic::Float64,
                                  kffi::Float64,
                                  ωad::Float64,
                                  kad::Float64)

            @vscontrol_ports

            new(kpv,
                kiv,
                kffv,
                rv,
                lv,
                kpc,
                kic,
                kffi,
                ωad,
                kad,
                6,
                [:ξ_d, :ξ_q, :γ_d, :γ_q, :ϕ_d, :ϕ_q],
                Ports(state_input, inner_input))

        end

end

"""
*  `kpc`::Float64 : current controller proportional gain
*  `kic`::Float64 : current controller integral gain
*  `kpv`::Float64 : voltage controller proportional gain
*  `kiv`::Float64 : voltage controller integral gain
*  `ω_ad`::Float64 : active damping filter cutoff frequency (rad/sec)
*  `kad`::Float64 : active damping gain
*  `lv`::Float64 : virtual inductance
*  `rv`::Float64 : virtual resistance
"""
mutable struct DirectControl <: VSControl
    val::Float64
    n_states::Int64
    states::Vector{Symbol}
    ports::Ports

        function DirectControl(val::Float64)

            @vscontrol_ports

            new(0,
                Vector{Symbol}(),
                Ports(state_input, inner_input))

        end

end
