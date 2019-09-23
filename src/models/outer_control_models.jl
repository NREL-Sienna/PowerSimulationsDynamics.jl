function mdl_outer_ode!(device_states,
                        output_ode,
                        f0,
                        device::DynInverter{C,VirtualInertiaQdroop{VirtualInertia,ReactivePowerDroop},VC,DC,P,F}) where {C <: Converter,
                                                                VC<: VSControl,
                                                                DC<: DCSource,
                                                                P <: FrequencyEstimator,
                                                                F <: Filter}


      #Obtain external states inputs for component
      external_ix = device.input_port_mapping[device.outercontrol]
      vpll_d = device_states[external_ix[1]]
      vpll_q = device_states[external_ix[2]]
      ϵ_pll = device_states[external_ix[3]]
      vod = device_states[external_ix[4]]
      voq = device_states[external_ix[5]]
      iod = device_states[external_ix[6]]
      ioq = device_states[external_ix[7]]

      #Obtain inner variables for component
      ω_pll =  device.inner_vars[ω_freq_estimator_var]

      #Get Active Power Controller parameters
      Ta = device.outercontrol.active_power.Ta #VSM Inertia constant
      kd = device.outercontrol.active_power.kd #VSM damping constant
      kω = device.outercontrol.active_power.kω #Frequency droop gain
      ωb = device.outercontrol.active_power.ωb #Rated angular frequency

      #Get Reactive Power Controller parameters
      kq = device.outercontrol.reactive_power.kq #Reactive power droop gain
      ωf = device.outercontrol.reactive_power.ωf #Reactive power filter cutoff frequency

      #Obtain external parameters
      kp_pll = device.freq_estimator.kp_pll
      ki_pll = device.freq_estimator.ki_pll
      p_ref = device.P_ref
      ω_ref = device.ω_ref
      V_ref = device.V_ref
      q_ref = device.Q_ref
      ωg = 1.0

      #Obtain indices for component w/r to device
      local_ix = device.local_state_ix[device.outercontrol]

      #Define internal states for frequency estimator
      internal_states = @view device_states[local_ix]
      δω_vsm = internal_states[1]
      δθ_vsm = internal_states[2]
      qm = internal_states[3]

      #Compute 3 states ODEs
      output_ode[local_ix[1]] = (- iod*vod/Ta
                                    - ioq*voq/Ta
                                    + kd*kp_pll*atan(vpll_q/vpll_d)/Ta
                                    + kd*ki_pll*ϵ_pll/Ta
                                    - (kd+kω)*δω_vsm/Ta
                                    + p_ref/Ta
                                    + kω*ω_ref/Ta
                                    - kω*ωg/Ta)
      output_ode[local_ix[2]] = ωb*δω_vsm
      output_ode[local_ix[3]] = (- ωf*ioq*vod
                                 + ωf*iod*voq
                                 - ωf*qm)

    #Update inner vars
    device.inner_vars[δdqRI_var] = δθ_vsm
    device.inner_vars[ω_control_var] = δω_vsm + 1.0
    device.inner_vars[v_control_var] = V_ref + kq*(q_ref - qm)
end
