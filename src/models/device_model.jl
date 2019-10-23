function device_model!(x,
                      output_ode::Vector{Float64},
                      voltage_r,
                      voltage_i,
                      current_r,
                      current_i,
                      ix_range::UnitRange{Int64},
                      ode_range::UnitRange{Int64},
                      controls,
                      device::DynG,
                      sys::DynamicSystem) where {DynG <: DynGenerator}
       #Obtain local device states
    n_states = total_device_states(device)
    device_states = @view x[ix_range]
       #Obtain references
    sys_Sbase = get_sys_base(sys)
    sys_f = get_sys_f(sys)

    #Update Voltage data
    device.inner_vars[VR_gen_var] = voltage_r[1]
    device.inner_vars[VI_gen_var] = voltage_i[1]

       #Obtain ODEs and Mechanical Power for Turbine Governor
    mdl_tg_ode!(device_states,
                view(output_ode, ode_range),
                device)
       #Obtain ODEs for AVR
    mdl_pss_ode!(device_states,
                 view(output_ode, ode_range),
                 device)
       #Obtain ODEs for AVR
    mdl_avr_ode!(device_states,
                view(output_ode, ode_range),
                device)
                   #Obtain ODEs for Machine
    mdl_machine_ode!(device_states,
                    view(output_ode, ode_range),
                    current_r,
                    current_i,
                    sys_Sbase,
                    sys_f,
                    device)
                       #Obtain ODEs for Shaft
    mdl_shaft_ode!(device_states,
                   view(output_ode, ode_range),
                   sys_f,
                   device)

    return

end


function device_model!(voltage_r,
                       voltage_i,
                       current_r,
                       current_i,
                       device::StaticSource,
                       sys::DynamicSystem)

    mdl_source!(voltage_r, voltage_i, current_r, current_i, device, sys)

    return
end

function device_model!(voltage_r,
                       voltage_i,
                       current_r,
                       current_i,
                       device::PSY.PowerLoad,
                       sys::DynamicSystem)

    mdl_Zload!(voltage_r, voltage_i, current_r, current_i, device, sys)

    return
end

function device_model!(x,
                        output_ode::Vector{Float64},
                        voltage_r,
                        voltage_i,
                        current_r,
                        current_i,
                        ix_range::UnitRange{Int64},
                        ode_range::UnitRange{Int64},
                        controls,
                        device::DynI,
                        sys::DynamicSystem) where {DynI <: DynInverter}

        #Obtain local device states
        n_states = total_device_states(device)
        device_states = @view x[ix_range]

        #Obtain references
        Sbase = get_sys_base(sys)
        sys_f = get_sys_f(sys)

        #Update Voltage data
        device.inner_vars[VR_inv_var] = voltage_r[1]
        device.inner_vars[VI_inv_var] = voltage_i[1]
        #device.inner_vars[Vh_var] = sqrt(voltage_r[1]^2 + voltage_i[1]^2)

        #Update V_ref
        V_ref = get_inverter_Vref(device)
        device.inner_vars[v_control_var] = V_ref

        #Obtain ODES for DC side
        mdl_DCside_ode!(device)

        #Obtain ODEs for OuterLoop
        mdl_outer_ode!(device_states,
                       view(output_ode, ode_range),
                       sys_f,
                       device)

#Obtain ODEs for PLL
mdl_freq_estimator_ode!(device_states,
                        view(output_ode, ode_range),
                        sys_f,
                        device)

#Obtain inner controller ODEs and modulation commands
mdl_VScontrol_ode!(device_states,
                   view(output_ode, ode_range),
                   device)

#Obtain converter relations
mdl_converter_ode!(device)

#Obtain ODEs for output filter
mdl_filter_ode!(device_states,
                        view(output_ode, ode_range),
                        current_r,
                        current_i,
                        Sbase,
                        sys_f,
                        device)




        return
end
