# This is where we do some type pyracy on the PSY types
get_n_buses(sys::PSY.System) = length(sys.bus_numbers)

function get_injectors_with_dynamics(sys::PSY.System)
    return PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing && PSY.get_available(x),
    )
end
