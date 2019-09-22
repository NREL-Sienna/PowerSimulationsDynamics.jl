abstract type Inverter end# <: DynInjection end

  """
  Parameters of the D'Arco Inverter Model

  # Constructor
  ```julia
  inverter_DAIM(Xd_p, f, H, D, MVABase, ω_ref)
  ```

  #Arguments
*  `p_ref`::Float64 : active power reference
*  `ω_ref`::Float64 : inverter frequency setpoint
*  `ω_ref`::Float64 : reference inverter frequency setpoint

*  `ωg``::Float64 : grid frequency

*  `q_ref`::Float64 : reactive power reference

*  `v_ref`::Float64 : output voltage reference

*  `lg`::Float64 : grid inductance
*  `rg`::Float64 : grid resistance
*  `vg`::Float64 : grid voltage

*  `kffv`::Float64 : Binary enable voltage feed-forward in current control
*  `kffi`::Float64 : Binary enable current feed-forwrad in current control

  """

mutable struct DAIM <: Inverter
    Ta::Float64
                  kd::Float64
                  kω::Float64
                  v_rated::Float64
                  s_rated::Float64
                  MVABase::Float64
                  ωb::Float64
                  ωg::Float64
                  kpc::Float64
                  kic::Float64
                  kpv::Float64
                  kiv::Float64
                  kq::Float64
                  ωf::Float64
                  lf::Float64
                  rf::Float64
                  cf::Float64
                  lg::Float64
                  rg::Float64
                  vg::Float64
                  ω_ad::Float64
                  kad::Float64
                  lv::Float64
                  rv::Float64
                  ω_lp::Float64
                  kp_pll::Float64
                  ki_pll::Float64
                  kffv::Float64
                  kffi::Float64
  n_states::Int64
  states::Vector{Symbol}

    function DAIM(Ta::Float64,
                  kd::Float64,
                  kω::Float64,
                  v_rated::Float64,
                  s_rated::Float64,
                  MVABase::Float64,
                  ωb::Float64,
                  ωg::Float64,
                  kpc::Float64,
                  kic::Float64,
                  kpv::Float64,
                  kiv::Float64,
                  kq::Float64,
                  ωf::Float64,
                  lf::Float64,
                  rf::Float64,
                  cf::Float64,
                  lg::Float64,
                  rg::Float64,
                  vg::Float64,
                  ω_ad::Float64,
                  kad::Float64,
                  lv::Float64,
                  rv::Float64,
                  ω_lp::Float64,
                  kp_pll::Float64,
                  ki_pll::Float64,
                  kffv::Float64,
                  kffi::Float64)

            n_states = 19
            states = [:δω_vsm, :δθ_vsm, :vod, :voq, :icvd, :icvq, :ξ_d, :ξ_q,
                     :γ_d, :γ_q, :iod, :ioq, :ϕ_d, :ϕ_q, :vpll_d, :vpll_q,
                     :ϵ_pll, :δθ_pll, :qm]
        new(Ta,
            kd,
            kω,
            v_rated,
            s_rated,
            MVABase,
            ωb,
            ωg,
            kpc,
            kic,
            kpv,
            kiv,
            kq,
            ωf,
            lf,
            rf,
            cf,
            lg,
            rg,
            vg,
            ω_ad,
            kad,
            lv,
            rv,
            ω_lp,
            kp_pll,
            ki_pll,
            kffv,
            kffi,
            n_states,
            states)
    end
end
