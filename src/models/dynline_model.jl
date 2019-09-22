function mdl_line_ode!(device_states,
                        output_ode,
                        V_r_from,
                        V_i_from,
                        V_r_to,
                        V_i_to,
                        current_r_from,
                        current_i_from,
                        current_r_to,
                        current_i_to,
                        dv_from,
                        dv_to,
                        sys_f,
                        branch)

    L = branch.x
    R = branch.r
    ω_b = sys_f*2*π
    c_from = branch.b.from
    c_to = branch.b.to


    current_r_from[1] -= (1.0/ω_b)*c_from*dv_from[1] - c_from*V_i_from
    current_i_from[1] -= (1.0/ω_b)*c_from*dv_from[2] + c_from*V_r_from
    current_r_to[1] -= (1.0/ω_b)*c_to*dv_to[1] - c_to*V_i_to
    current_i_to[1] -= (1.0/ω_b)*c_to*dv_to[2] + c_to*V_r_to


    Il_r = device_states[1]
    Il_i = device_states[2]
    output_ode[1] = (ω_b/L)*((V_r_from - V_r_to) - (R*Il_r - L*Il_i))
    output_ode[2] = (ω_b/L)*((V_i_from - V_i_to) - (R*Il_i + L*Il_r))

    current_r_from[1] -= Il_r
    current_i_from[1] -= Il_i
    current_r_to[1] += Il_r
    current_i_to[1] += Il_i

end
