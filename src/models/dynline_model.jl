mutable struct DynamicLine{T<:PSY.ACBranch} <: PSY.Device
    branch_data::T
    n_states::Int64
    states

    function DynamicLine(branch::T) where T <: PSY.ACBranch
        n_states = 2
        states = [
            :Il_R
            :Il_I
        ]
        new{T}(branch, n_states, states)
    end

end

"""Get Line name."""
IS.get_name(value::DynamicLine) = value.branch_data.name
"""Get Line available."""
PSY.get_available(value::DynamicLine) = value.branch_data.available
"""Get Line activepower_flow."""
PSY.get_activepower_flow(value::DynamicLine) = value.branch_data.activepower_flow
"""Get Line reactivepower_flow."""
PSY.get_reactivepower_flow(value::DynamicLine) = value.branch_data.reactivepower_flow
"""Get Line arc."""
PSY.get_arc(value::DynamicLine) = value.branch_data.arc
"""Get Line r."""
PSY.get_r(value::DynamicLine) = value.branch_data.r
"""Get Line x."""
PSY.get_x(value::DynamicLine) = value.branch_data.x
"""Get Line b."""
PSY.get_b(value::DynamicLine) = value.branch_data.b
"""Get Line rate."""
PSY.get_rate(value::DynamicLine) = value.branch_data.rate
"""Get Line anglelimits."""
PSY.get_anglelimits(value::DynamicLine) = value.branch_data.anglelimits
"""Get Line services."""
PSY.get_services(value::DynamicLine) = value.branch_data.services
"""Get Line ext."""
PSY.get_ext(value::DynamicLine) = value.branch_data.ext
"""Get Line _forecasts."""
PSY.get__forecasts(value::DynamicLine) = value.branch_data._forecasts
"""Get Line internal."""
PSY.get_internal(value::DynamicLine) = value.branch_data.internal
PSY.get_states(value::DynamicLine) = value.states
PSY.get_n_states(value::DynamicLine) = value.n_states

function make_dynamic_branch!(branch::T, sys::PSY.System) where T <: PSY.ACBranch
    PSY.remove_component!(sys, branch)
    PSY.add_component!(sys, DynamicLine(branch))
    return
end

function mdl_line_ode!(
    device_states,
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
    branch,
)

    L = branch.x
    R = branch.r
    ω_b = sys_f * 2 * π
    c_from = branch.b.from
    c_to = branch.b.to

    current_r_from[1] -= (1.0 / ω_b) * c_from * dv_from[1] - c_from * V_i_from[1]
    current_i_from[1] -= (1.0 / ω_b) * c_from * dv_from[2] + c_from * V_r_from[1]
    current_r_to[1] -= (1.0 / ω_b) * c_to * dv_to[1] - c_to * V_i_to[1]
    current_i_to[1] -= (1.0 / ω_b) * c_to * dv_to[2] + c_to * V_r_to[1]

    Il_r = device_states[1]
    Il_i = device_states[2]
    output_ode[1] = (ω_b / L) * ((V_r_from[1] - V_r_to[1]) - (R * Il_r - L * Il_i))
    output_ode[2] = (ω_b / L) * ((V_i_from[1] - V_i_to[1]) - (R * Il_i + L * Il_r))

    current_r_from[1] -= Il_r
    current_i_from[1] -= Il_i
    current_r_to[1] += Il_r
    current_i_to[1] += Il_i

end
