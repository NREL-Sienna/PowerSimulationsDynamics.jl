function mdl_DAIM_ode(internal_states,
                      external_states,
                      inputs::Vector{Float64},
                      outputs,
                      p::Inverter)#,
                      #gen_ix::Dict{Symbol,Int64})
    #Define inner states for device
    δω_vsm = internal_states[1]
    δθ_vsm = internal_states[2]
       vod = internal_states[3]
       voq = internal_states[4]
      icvd = internal_states[5]
      icvq = internal_states[6]
       ξ_d = internal_states[7]
       ξ_q = internal_states[8]
       γ_d = internal_states[9]
       γ_q = internal_states[10]
       iod = internal_states[11]
       ioq = internal_states[12]
       ϕ_d = internal_states[13]
       ϕ_q = internal_states[14]
    vpll_d = internal_states[15]
    vpll_q = internal_states[16]
     ϵ_pll = internal_states[17]
    δθ_pll = internal_states[18]
        qm = internal_states[19]

    #Define inputs for device
    p_ref = controls[1]
    q_ref = controls[2]
    ω_ref = controls[3]
    v_ref = controls[4]
    ω_grd = controls[5]

    #Define external states for device
    V_R = external_states[1]
    V_I = external_states[2]

    V_dq = ri_dq(δθ_vsm)*[V_R; V_I]

    ### function([internal_states],
    ###          [external_states],
    ###          [inputs],[outputs],params)
    inverter_ODE = vcat(
        filter_LCL([icvd,icvq,iod,ioq,vod,voq],
                   [vcvd,vcvq,vgd,vgq],
                   [],[],p),
        SRF_current_ctrl_extDC([ϕ_d,ϕ_q,γ_d,γ_q],
                   [vod,voq,icvd,icvq,vdc],
                   [ω_vsm,icvd_ref,icvq_ref],[md,mq],p),
        SRF_voltage_VZ_ctrl([ξ_d,ξ_q],
                   [vod,voq,iod,ioq],
                   [ω_vsm,v_vrefr],[icvd_ref,icvq_ref],p),
        reactive_power_droop_ctrl([qm],
                   [vod,voq,iod,ioq],
                   [v_ref,q_ref],[v_refr],p),
        active_power_VSM_droop_ctrl([dδω_vsm,dδθ_vsm],
                   [vod,voq,iod,ioq],
                   [p_ref,ω_ref,ω_pll],[ω_vsm],p),
        PLL_dq_PI_ctrl([vpll_d,vpll_q,ϵ_pll,δθ_pll],
                   [vod,voq,δθ_vsm],
                   [ω_grd],[ω_pll],p)
        )
    inverter_DAE = inverter_avg_extDC([vcvd,vcvq,idc],
                   [icvd,icvq,vdc],
                   [md,mq],[],p)

    I_RI_vec = p.MVABase*dq_ri(δθ_vsm)*[iod; ioq]

    return inverter_ODE, [I_RI_vec,inverter_DAE]
end

function filter_LCL(internal_states,
                    external_states,
                    inputs,
                    outputs,
                    params::DAIM)

    #Define internal states for filter
    icvd = internal_states[1]
    icvq = internal_states[2]
     iod = internal_states[3]
     ioq = internal_states[4]
     vod = internal_states[5]
     voq = internal_states[6]

    #Define external states for filter
    vcvd = external_states[1]
    vcvq = external_states[2]
     vgd = external_states[3]
     vgq = external_states[4]

    #Inputs (control signals)
    wg = inputs[1]

    #Get parameters
    ωb = params.ωb
    lf = params.lf
    rf = params.rf
    cf = params.cf
    lg = params.lg

    #ODEs
    #Inverter Output Inductor (internal state)
    dicvd_dt = ( ωb/lf*vcvd
               - ωb/lf*vod
               - ωb*rf/lf*icvd
               - ωb*ωg*icvq )
    dicvq_dt = ( ωb/lf*vcvq
               - ωb/lf*voq
               - ωb*rf/lf*icvq
               + ωb*ωg*icvq )
    #LCL Capacitor (internal state)
    dvod_dt = ( ωb/cf*icvd
              - ωb/cf*iod  #i_gd was specified; use equivalent i_od
              - ωb*ωg*voq )
    dvoq_dt = ( ωb/cf*icvq
              - ωb/cf*ioq  #i_gq was specified; use equivalent i_oq
              + ωb*ωg*vod )
    #Grid Inductance (internal state)
    diod_dt = ( ωb*vod/lg
              - ωb*rg*iod/lg
              + ωb*ωg*ioq
              - ωb*vgd/lg)
    dioq_dt = ( ωb*voq/lg
              - ωb*ωg*iod
              - ωb*rg*ioq/lg
              + ωb*vgq/lg)

    #Return ODE
    return [dicvd_dt, dicvq_dt, dvod_dt, dvoq_dt, diod_dt, dioq_dt]
end
function inverter_avg_fixedDC(internal_states,
                              external_states,
                              inputs,
                              outputs,
                              params::DAIM)
    #ALL ALGEBRAIC
    #Get internal states for avg inverter model
    vcvd = internal_states[1]
    vcvq = internal_states[2]

    #Get parameters
    vdc = 1.0

    #Compute algebraic
    vcvd = md*vdc
    vcvq = mq*vdc

    #Outputs (control signals)
    outputs[1] = md
    outputs[2] = mq
    #Return DAE
    return [vcvd,vcvq]
end
function inverter_avg_extDC(internal_states,
                            external_states,
                            inputs,
                            outputs,
                            params::DAIM)
    #ALL ALGEBRAIC
    #Get internal states for avg inverter model
    vcvd = internal_states[1]
    vcvq = internal_states[2]
     idc = internal_states[3]

    #Get external states for avg inverter model
    icvd = external_states[1]
    icvq = external_states[2]
     vdc = external_states[3]

    #Inputs (control signals)
    md = inputs[1]
    mq = inputs[2]

    #Compute algebraic
    vcvd = md*vdc
    vcvq = mq*vdc
     idc = (vcvd*icvd+vcvq*icvq)/vdc

    #Return DAE
    return [vcvd,vcvq,idc]
end
function SRF_current_ctrl_extDC(internal_states,
                                external_states,
                                inputs,
                                outputs,
                                params::DAIM)
    #Get internal states for current controller
    ϕ_d = internal_states[1]
    ϕ_q = internal_states[2]
    γ_d = internal_states[3]
    γ_q = internal_states[4]

    #Get external states for current controller
     vod = external_states[1]
     voq = external_states[2]
    icvd = external_states[3]
    icvd = external_states[4]
     vdc = external_states[5]

    #Inputs (control signals)
       ω_vsm = inputs[1]
    icvd_ref = inputs[2]
    icvq_ref = inputs[3]

    #Get parameters
     kpc = params.kpc
     kic = params.kic
    kffv = params.kffv
      lf = params.lf
     ωad = params.ωad
     kad = params.kad

    #Computation but not DAE
    vad_d = kad*(vod-ϕ_d)
    vad_q = kad*(voq-ϕ_q)

    vcvd_ref = ( kpc*(icvd_ref-icvd)
               + kic*γ_d
               + ω_vsm*lf*icvq #j product - change axis
               + kffv*vod
               - vad_d )
    vcvq_ref = ( kpc*(icvq_ref-icvq)
               + kic*γ_q
               - ω_vsm*lf*icvd #j product - change axis
               + kffv*voq
               - vad_q )
    #Output Control Signal
    md = vcvd_ref/vdc
    mq = vcvq_ref/vdc

    #ODEs
    #Active Damping LPF (internal state)
    dϕ_d_dt = ωad*vod - ωad*ϕ_d
    dϕ_q_dt = ωad*voq - ωad*ϕ_q
    #PI Integrator (internal state)
    dγ_d_dt = icvd_ref - icvd
    dγ_q_dt = icvq_ref - icvq

    #Outputs (control signals)
    outputs[1] = md
    outputs[2] = mq
    #Return ODE
    return [dϕ_d_dt, dϕ_q_dt, dγ_d_dt, dγ_q_dt]
end
function SRF_voltage_VZ_ctrl(internal_states,
                                   external_states,
                                   inputs,
                                   outputs,
                                   params::DAIM)
    #Get internal states for voltage controller
    ξ_d = internal_states[1]
    ξ_q = internal_states[2]

    #Get external states for voltage controller
    vod = external_states[1]
    voq = external_states[2]
    iod = external_states[3]
    ioq = external_states[4]

    #Inputs (control signals)
     ω_vsm = inputs[1]
    v_refr = inputs[2]

    #Get parameters
     kpv = params.kpv
     kiv = params.kiv
    kffi = params.kffi
      cf = params.cf
      rv = params.rv
      lv = params.lv

    #Virtual Impedance - Computation but not DAE
    vod_ref = ( v_refr - rv*iod + ω_vsm*lv*ioq )
    voq_ref = (          rv*ioq - ω_vsm*lv*iod )
    #Output Control Signal
    icvd_ref = ( kpv*(vod_ref-vod)
               + kiv*ξ_d
               + cf*ω_vsm*voq #j product - change axis
               + kffi*iod )

    icvq_ref = ( kpv*(voq_ref-voq)
               - kiv*ξ_q
               - cf*ω_vsm*vod #j product - change axis
               + kffi*ioq )

    #ODE
    #PI Integrator (internal state)
    dξ_d_dt = (vod_ref-vod)
    dξ_q_dt = (voq_ref-voq)

    #Outputs (control signals)
    outputs[1] = icvd_ref
    outputs[2] = icvq_ref
    #Return ODE
    return [dξ_d_dt, dξ_q_dt]
end
function reactive_power_droop_ctrl(internal_states,
                                external_states,
                                inputs,
                                outputs,
                                params::DAIM)
    #Get internal states for reactive power
    qm = internal_states[1]

    #Get external states for reactive power
    vod = external_states[1]
    voq = external_states[2]
    iod = external_states[3]
    ioq = external_states[4]

    #Inputs (control signals)
    v_ref = inputs[1]
    q_ref = inputs[2]

    #Get parameters
    ωf = params.ωf
    kq = params.kq

    #Output Control Signal
    v_refr = v_ref - kq*(q_ref-qm)

    #ODE
    #Reactive Power LPF (internal state)
    dqm_dt = -ωf*qm+ωf*(vod*ioq+voq*iod)

    #Outputs (control signals)
    outputs[1] = v_refr
    #Return ODE
    return [dqm_dt]
end
function active_power_VSM_droop_ctrl(internal_states,
                                external_states,
                                inputs,
                                outputs,
                                params::DAIM)
    #Get internal states for active power
    dδω_vsm = internal_states[1]
    dδθ_vsm = internal_states[2]

    #Get external states for active power
    vod = external_states[1]
    voq = external_states[2]
    iod = external_states[3]
    ioq = external_states[4]

    #Inputs (control signals)
    p_ref = inputs[1]
    ω_ref = inputs[2]
    ω_pll = inputs[3] # params.ωg

    #Get parameters
    Ta = params.Ta
    kd = params.kd
    kω = params.kω
    ωb = params.ωb

    #Compute Outputs
    ω_vsm = δω_vsm + w_grd

    #ODE
    #VSM Frequency Deviation (internal state)
    dδω_vsm_dt = ( p_ref/Ta
                 - (iod*vod+ioq*voq)/Ta
                 - kd*(ω_vsm-ω_pll)/Ta
                 - kω*(ω_vsm-ω_ref)/Ta)
    #VSM Angle Displacement (internal state)
    dδθ_vsm_dt = ωb*δω_vsm

    #Outputs (control signals)
    outputs[1] = ω_vsm
    #Return ODE
    return [dδω_vsm_dt, dδθ_vsm_dt]
end
function PLL_dq_PI_ctrl(internal_states,
                        external_states,
                        inputs,
                        outputs,
                        params::DAIM)

    #Get internal states for PLL
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
    ϵ_pll = internal_states[3]
    δθ_pll = internal_states[4]

    #Get external states for PLL
    vod = external_states[1]
    voq = external_states[2]
    δθ_vsm = external_states[3]

    #Inputs (control signals)
    ω_grd = controls[1] #params.ωg

    #Get parameters
    ω_lp = params.ω_lp
    kp_pll = params.kp_pll
    ki_pll = params.ki_pll
    ωb = params.ωb

    #ODE
    #Output Voltage LPF (internal state)
    dvpll_d_dt =  (ω_lp*vod*cos(δθ_pll-δθ_vsm)
                 + ω_lp*voq*sin(δθ_pll-δθ_vsm)
                 - ω_lp*vpll_d)
    dvpll_q_dt = (- ω_lp*vod*sin(δθ_pll-δθ_vsm)
                  + ω_lp*voq*cos(δθ_pll-δθ_vsm)
                  - ω_lp*vpll_q)
    #PI Integrator (internal state)
    dϵ_pll_dt = atan(vpll_q/vpll_d)
    #PLL Frequency Deviation (internal state)
    dδθ_pll =  ( ωb*kp_pll*atan(vpll_q/vpll_d)
               + ωb*ki_pll*ϵ_pll)

    #Outputs (control signals)
    outputs[1] = ω_pll
    #Return ODEs
    return [dvpll_d_dt, dvpll_q_dt, dϵ_pll_dt, dδθ_pll]
end
