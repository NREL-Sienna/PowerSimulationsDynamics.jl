function mdl_DAIM_ode(internal_states,
                      external_states,
                      controls::Vector{Float64},
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

    inverter_ODE = vcat(vsm_inertia(δω_vsm, [iod, vod, ioq, voq, vpll_d, vpll_q, ϵ_pll], [p_ref, ω_ref, ω_grd], p),
                        voltage_control([vod, voq, icvd, icvq, ξ_d, ξ_q], [iod, ioq, ϕ_d, ϕ_q, γ_d, γ_q, qm, δω_vsm], [q_ref, v_ref, ω_grd], p),
                        current_control_extref([iod, ioq, ϕ_d, ϕ_q], [vod, voq, icvd, icvq, ξ_d, ξ_q, qm, δθ_vsm, δω_vsm, V_dq[2], V_dq[1]], [q_ref, v_ref, ω_grd], p),
                        PLL([vpll_d, vpll_q, ϵ_pll, δθ_pll], [vod, voq, δθ_vsm], ω_grd, p),
                        reactive_power_droop(qm, [iod, vod, ioq, voq], p))

    I_RI_vec = p.MVABase*dq_ri(δθ_vsm)*[iod; ioq]

    return inverter_ODE, I_RI_vec
end

function vsm_inertia(internal_states,
                     external_states,
                     controls,
                     params::DAIM)

    #Define inner states for inertia
    δω_vsm = internal_states[1]

    #Define external states for inerita
    iod = external_states[1]
    vod = external_states[2]
    ioq = external_states[3]
    voq = external_states[4]
    vpll_d = external_states[5]
    vpll_q = external_states[6]
    ϵ_pll = external_states[7]

    #Get parameters
    Ta = params.Ta
    kd = params.kd
    kp_pll = params.kp_pll
    ki_pll = params.ki_pll
    kω = params.kω
    p_ref = controls[1]
    ω_ref = controls[2]
    ωg = controls[3] # params.ωg
    ωb = params.ωb

    #Compute ODEs
    dδω_vsm_dt = (- iod*vod/Ta
                 - ioq*voq/Ta
                 + kd*kp_pll*atan(vpll_q/vpll_d)/Ta
                 + kd*ki_pll*ϵ_pll/Ta
                 - (kd+kω)*δω_vsm/Ta
                 + p_ref/Ta
                 + kω*ω_ref/Ta
                 - kω*ωg/Ta)

    dδθ_vsm_dt = ωb*δω_vsm

    #Return ODEs
    return [dδω_vsm_dt, dδθ_vsm_dt]
end

function voltage_control(internal_states,
                         external_states,
                         controls,
                         params::DAIM)

    #Define inner states for inertia
    vod = internal_states[1]
    voq = internal_states[2]
    icvd = internal_states[3]
    icvq = internal_states[4]
    ξ_d = internal_states[5]
    ξ_q = internal_states[6]

    #Define external states for inerita
    iod = external_states[1]
    ioq = external_states[2]
    ϕ_d = external_states[3]
    ϕ_q = external_states[4]
    γ_d = external_states[5]
    γ_q = external_states[6]
    qm = external_states[7]
    δω_vsm = external_states[8]

    #Get parameters
    kpv = params.kpv
    kiv = params.kiv
    kq = params.kq
    ωb = params.ωb
    #ωg = params.ωg
    lf = params.lf
    rf = params.rf
    cf = params.cf
    kad = params.kad
    kffv = params.kffv
    kffi = params.kffi
    kpc = params.kpc
    kic = params.kic
    rv = params.rv
    lv = params.lv
    q_ref = controls[1]
    v_ref = controls[2]
    ωg = controls[3]

    #Compute ODEs
    dvod_dt =   (ωb*ωg*voq
              + ωb*icvd/cf
              - ωb*iod/cf)

    dvoq_dt = (- ωb*ωg*vod
              + ωb*icvq/cf
              - ωb*ioq/cf)

    dicvd_dt = (  (ωb*(kffv-1.0-kad-kpc*kpv))*vod/lf
               - ωb*cf*kpc*ωg*voq/lf
               - ωb*(kpc+rf)*icvd/lf
               + ωb*kic*γ_d/lf
               + ωb*kpc*(kffi-kpv*rv)*iod/lf
               + ωb*kpc*kpv*lv*ωg*ioq/lf
               + ωb*kad*ϕ_d/lf
               + ωb*kiv*kpc*ξ_d/lf
               - ωb*kpc*kpv*kq*qm/lf
               - ωb*icvq*δω_vsm
               + ωb*kpc*kpv*lv*ioq*δω_vsm/lf
               - ωb*cf*kpc*voq*δω_vsm/lf
               + ωb*kpc*kpv*kq*q_ref/lf
               + ωb*kpc*kpv*v_ref/lf )

    dicvq_dt = ( ωb*cf*kpc*ωg*vod/lf
               + (ωb*(kffv-1.0-kad-kpc*kpv))*voq/lf
               - ωb*(kpc+rf)*icvq/lf
               + ωb*kic*γ_q/lf
               - ωb*kpc*kpv*lv*ωg*iod/lf
               + ωb*kpc*(kffi-kpv*rv)*ioq/lf
               + ωb*kad*ϕ_q/lf
               + ωb*kiv*kpc*ξ_q/lf
               + ωb*icvd*δω_vsm
               - ωb*kpc*kpv*lv*iod*δω_vsm/lf
               + ωb*cf*kpc*vod*δω_vsm/lf )

    dξ_d_dt = (- vod
              - rv*iod
              + lv*ωg*ioq
              - kq*qm
              + lv*ioq*δω_vsm
              + kq*q_ref
              + v_ref )

    dξ_q_dt = (- voq
              - lv*ωg*iod
              - rv*ioq
              - lv*iod*δω_vsm )

    #Return ODEs
    return [dvod_dt, dvoq_dt, dicvd_dt, dicvq_dt, dξ_d_dt, dξ_q_dt]
end

function current_control_extref(internal_states,
                                external_states,
                                controls,
                                params::DAIM)

    #Define inner states for inertia
    iod = internal_states[1]
    ioq = internal_states[2]
    ϕ_d = internal_states[3]
    ϕ_q = internal_states[4]

    #Define external states for inerita
    vod = external_states[1]
    voq = external_states[2]
    icvd = external_states[3]
    icvq = external_states[4]
    ξ_d = external_states[5]
    ξ_q = external_states[6]
    qm = external_states[7]
    δθ_vsm = external_states[8]
    δω_vsm = external_states[9]
    vgd = external_states[10]
    vgq = external_states[11]

    #Get parameters
    kpv = params.kpv
    kiv = params.kiv
    ωb = params.ωb
    #ωg = params.ωg
    cf = params.cf
    kffi = params.kffi
    lg = params.lg
    rg = params.rg
    ω_ad = params.ω_ad
    rv = params.rv
    lv = params.lv
    kq = params.kq
    q_ref = controls[1]
    v_ref = controls[2]
    ωg = controls[3]
    vg = params.vg

    #Compute ODEs
    dγ_d_dt = (- kpv*vod
              - cf*ωg*voq
              - icvd
              + (kffi-kpv*rv)*iod
              + kpv*lv*ωg*ioq
              + kiv*ξ_d
              - kpv*kq*qm
              + kpv*lv*ioq*δω_vsm
              - cf*voq*δω_vsm
              + kpv*kq*q_ref
              + kpv*v_ref )

    dγ_q_dt =  ( cf*ωg*vod
             - kpv*voq
             - icvq
             - kpv*lv*ωg*iod
             + (kffi-kpv*rv)*ioq
             + kiv*ξ_q
             - kpv*lv*iod*δω_vsm
             + cf*vod*δω_vsm )

    diod_dt = ( ωb*vod/lg
             - ωb*rg*iod/lg
             + ωb*ωg*ioq
             - ωb*vgd/lg ) #EPSR 122 Paper eqn error: should be (-) term

    dioq_dt = (ωb*voq/lg
             - ωb*ωg*iod
             - ωb*rg*ioq/lg
             + ωb*vgq/lg)

    dϕ_d_dt = ω_ad*vod - ω_ad*ϕ_d

    dϕ_q_dt = ω_ad*voq - ω_ad*ϕ_q

    #Return ODEs
    return [dγ_d_dt, dγ_q_dt, diod_dt, dioq_dt, dϕ_d_dt, dϕ_q_dt]
end

function PLL(internal_states,
             external_states,
             controls,
             params::DAIM)

    #Define inner states for inertia
    vpll_d = internal_states[1]
    vpll_q = internal_states[2]
    ϵ_pll = internal_states[3]
    δθ_pll = internal_states[4]

    #Define external states for inerita
    vod = external_states[1]
    voq = external_states[2]
    δθ_vsm = external_states[3]

    #Get parameters
    ω_lp = params.ω_lp
    kp_pll = params.kp_pll
    ki_pll = params.ki_pll
    ωb = params.ωb
    ωg = controls[1] #params.ωg

    #Compute ODEs
    dvpll_d_dt =  (ω_lp*vod*cos(δθ_pll-δθ_vsm)
                 + ω_lp*voq*sin(δθ_pll-δθ_vsm)
                 - ω_lp*vpll_d)

    dvpll_q_dt = (- ω_lp*vod*sin(δθ_pll-δθ_vsm)
                  + ω_lp*voq*cos(δθ_pll-δθ_vsm)
                  - ω_lp*vpll_q)

    dϵ_pll_dt = atan(vpll_q/vpll_d)

    dδθ_pll =  (ωb*kp_pll*atan(vpll_q/vpll_d)
            + ωb*ki_pll*ϵ_pll)

    #Return ODEs
    return [dvpll_d_dt, dvpll_q_dt, dϵ_pll_dt, dδθ_pll]
end

function reactive_power_droop(internal_states,
                              external_states,
                              params::DAIM)

    #Define inner states for inertia
    qm = internal_states[1]

    #Define external states for inerita
    iod = external_states[1]
    vod = external_states[2]
    ioq = external_states[3]
    voq = external_states[4]

    #Get parameters
    ωf = params.ωf

    #Compute ODEs
    dqm_dt = (- ωf*ioq*vod
             + ωf*iod*voq
             - ωf*qm)

    #Return ODEs
    return [dqm_dt]
end
