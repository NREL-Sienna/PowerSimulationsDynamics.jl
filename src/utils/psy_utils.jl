# This is where we do some type pyracy on the PSY types
get_n_buses(sys::PSY.System) = length(sys.bus_numbers)

function get_injectors_with_dynamics(sys::PSY.System)
    return PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing && PSY.get_available(x),
    )
end

function get_dynamic_branches(sys::PSY.System)
    return PSY.get_components(PSY.DynamicBranch, sys, x -> PSY.get_available(x))
end

function _transform_all_lines!(sys::PSY.System)
    for br in PSY.get_components(PSY.DynamicBranch, sys)
        dyn_br = DynamicBranch(br)
        @debug "Converted $(PSY.get_name(dyn_br)) to DynamicBranch"
        add_component!(sys, dyn_br)
    end
end
