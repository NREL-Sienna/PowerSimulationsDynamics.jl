# This is where we do some type piracy on the PSY types
get_n_buses(sys::PSY.System) = length(sys.bus_numbers)

function _filter_function(x::T) where {T <: PSY.StaticInjection}
    if PSY.get_dynamic_injector(x) === nothing
        return false
    end

    if hasfield(T, :status)
        return PSY.get_status(x) * PSY.get_available(x)
    else
        return PSY.get_available(x)
    end
    error("Filtering of $(PSY.get_name(x)) failed")
end

function get_injectors_with_dynamics(sys::PSY.System)
    return PSY.get_components(x -> _filter_function(x), PSY.StaticInjection, sys)
end

function get_injection_without_dynamics(sys::PSY.System)
    return PSY.get_components(
        x ->
            PSY.get_dynamic_injector(x) === nothing &&
                PSY.get_available(x) &&
                !isa(x, PSY.ElectricLoad),
        PSY.StaticInjection,
        sys,
    )
end

function get_dynamic_branches(sys::PSY.System)
    return PSY.get_components(x -> PSY.get_available(x), PSY.DynamicBranch, sys)
end

function _transform_all_lines!(sys::PSY.System)
    for br in PSY.get_components(PSY.DynamicBranch, sys)
        dyn_br = DynamicBranch(br)
        @debug "Converted $(PSY.get_name(dyn_br)) to DynamicBranch"
        add_component!(sys, dyn_br)
    end
end

function transform_ybus_to_rectangular(
    ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int},
)
    # TODO: Improve performance here
    return hcat(vcat(real(ybus), -imag(ybus)), vcat(imag(ybus), real(ybus)))
end

function transform_branches_to_dynamic(sys::PSY.System, ::Type{T}) where {T <: PSY.ACBranch}
    for b in PSY.get_components(T, sys)
        dyn_branch = PSY.DynamicBranch(b)
        PSY.add_component!(sys, dyn_branch)
    end
    return
end

function _compute_total_load_parameters(load::PSY.StandardLoad)
    # Constant Power Data
    constant_active_power = PSY.get_constant_active_power(load)
    constant_reactive_power = PSY.get_constant_reactive_power(load)
    max_constant_active_power = PSY.get_max_constant_active_power(load)
    max_constant_reactive_power = PSY.get_max_constant_reactive_power(load)
    # Constant Current Data
    current_active_power = PSY.get_current_active_power(load)
    current_reactive_power = PSY.get_current_reactive_power(load)
    max_current_active_power = PSY.get_max_current_active_power(load)
    max_current_reactive_power = PSY.get_max_current_reactive_power(load)
    # Constant Admittance Data
    impedance_active_power = PSY.get_impedance_active_power(load)
    impedance_reactive_power = PSY.get_impedance_reactive_power(load)
    max_impedance_active_power = PSY.get_max_impedance_active_power(load)
    max_impedance_reactive_power = PSY.get_max_impedance_reactive_power(load)
    # Total Load Calculations
    active_power = constant_active_power + current_active_power + impedance_active_power
    reactive_power =
        constant_reactive_power + current_reactive_power + impedance_reactive_power
    max_active_power =
        max_constant_active_power + max_current_active_power + max_impedance_active_power
    max_reactive_power =
        max_constant_reactive_power +
        max_current_reactive_power +
        max_impedance_reactive_power
    return active_power, reactive_power, max_active_power, max_reactive_power
end

function transform_load_to_constant_impedance(load::PSY.StandardLoad)
    # Total Load Calculations
    active_power, reactive_power, max_active_power, max_reactive_power =
        _compute_total_load_parameters(load)
    # Set Impedance Power
    PSY.set_impedance_active_power!(load, active_power)
    PSY.set_impedance_reactive_power!(load, reactive_power)
    PSY.set_max_impedance_active_power!(load, max_active_power)
    PSY.set_max_impedance_reactive_power!(load, max_reactive_power)
    # Set everything else to zero
    PSY.set_constant_active_power!(load, 0.0)
    PSY.set_constant_reactive_power!(load, 0.0)
    PSY.set_max_constant_active_power!(load, 0.0)
    PSY.set_max_constant_reactive_power!(load, 0.0)
    PSY.set_current_active_power!(load, 0.0)
    PSY.set_current_reactive_power!(load, 0.0)
    PSY.set_max_current_active_power!(load, 0.0)
    PSY.set_max_current_reactive_power!(load, 0.0)
    return
end

function transform_load_to_constant_current(load::PSY.StandardLoad)
    # Total Load Calculations
    active_power, reactive_power, max_active_power, max_reactive_power =
        _compute_total_load_parameters(load)
    # Set Impedance Power
    PSY.set_current_active_power!(load, active_power)
    PSY.set_current_reactive_power!(load, reactive_power)
    PSY.set_max_current_active_power!(load, max_active_power)
    PSY.set_max_current_reactive_power!(load, max_reactive_power)
    # Set everything else to zero
    PSY.set_constant_active_power!(load, 0.0)
    PSY.set_constant_reactive_power!(load, 0.0)
    PSY.set_max_constant_active_power!(load, 0.0)
    PSY.set_max_constant_reactive_power!(load, 0.0)
    PSY.set_impedance_active_power!(load, 0.0)
    PSY.set_impedance_reactive_power!(load, 0.0)
    PSY.set_max_impedance_active_power!(load, 0.0)
    PSY.set_max_impedance_reactive_power!(load, 0.0)
    return
end

function transform_load_to_constant_power(load::PSY.StandardLoad)
    # Total Load Calculations
    active_power, reactive_power, max_active_power, max_reactive_power =
        _compute_total_load_parameters(load)
    # Set Impedance Power
    PSY.set_constant_active_power!(load, active_power)
    PSY.set_constant_reactive_power!(load, reactive_power)
    PSY.set_max_constant_active_power!(load, max_active_power)
    PSY.set_max_constant_reactive_power!(load, max_reactive_power)
    # Set everything else to zero
    PSY.set_current_active_power!(load, 0.0)
    PSY.set_current_reactive_power!(load, 0.0)
    PSY.set_max_current_active_power!(load, 0.0)
    PSY.set_max_current_reactive_power!(load, 0.0)
    PSY.set_impedance_active_power!(load, 0.0)
    PSY.set_impedance_reactive_power!(load, 0.0)
    PSY.set_max_impedance_active_power!(load, 0.0)
    PSY.set_max_impedance_reactive_power!(load, 0.0)
    return
end
