function initialize_outer!(device_states,
    device::PSY.DynamicInverter{
        C,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
        IC,
        DC,
        P,
        F,
    },
) where {
    C <: PSY.Converter,
    IC <: PSY.InnerControl,
    DC <: PSY.DCSource,
    P <: PSY.FrequencyEstimator,
    F <: PSY.Filter,
}

    #Obtain external states inputs for component
    external_ix = get_input_port_ix(
        device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    Vr_filter = device_states[external_ix[4]]
    Vi_filter = device_states[external_ix[5]]
    Ir_filter = device_states[external_ix[6]]
    Ii_filter = device_states[external_ix[7]]

    Vr_cnv = get_inner_vars(device)[Vr_cnv_var]
    Vi_cnv = get_inner_vars(device)[Vi_cnv_var]
    θ0_oc = angle(Vr_cnv + 1im*Vi_cnv)

    #Obtain additional expressions
    p_elec_out = Ir_filter * Vr_filter + Ii_filter * Vi_filter
    q_elec_out = -Ii_filter * Vr_filter + Ir_filter * Vi_filter

    #Update inner_vars
    get_inner_vars(device)[P_ES_var] = p_elec_out
    #Update states
    outer_ix = get_local_state_ix(
        device,
        PSY.OuterControl{PSY.VirtualInertia, PSY.ReactivePowerDroop},
    )
    outer_states = @view device_states[outer_ix]
    outer_states[1] = PSY.get_ω_ref(device) #ω
    outer_states[2] = θ0_oc #θ_oc
    outer_states[3] = q_elec_out #qm

    #Update inner vars
    get_inner_vars(device)[θ_oc_var] = θ0_oc
    get_inner_vars(device)[ω_oc_var] = PSY.get_ω_ref(device)
    #Update Q_ref. Initialization assumes q_ref = q_elec_out of PF solution
    PSY.get_ext(device)[CONTROL_REFS][P_ref_index] = p_elec_out
    PSY.set_P_ref!(PSY.get_active_power(PSY.get_outer_control(device)), p_elec_out)
    PSY.get_ext(device)[CONTROL_REFS][Q_ref_index] = q_elec_out
end