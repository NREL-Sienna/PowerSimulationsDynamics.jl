function power_flow_solution!(
    initial_guess::Vector{Float64},
    sys::PSY.System,
    inputs::SimulationInputs,
)
    res = PF.solve_ac_powerflow!(sys)
    if !res
        CRC.@ignore_derivatives @error("PowerFlow failed to solve")
        return BUILD_FAILED
    end
    bus_size = length(PSY.get_bus_numbers(sys))
    CRC.@ignore_derivatives @debug "Updating bus voltage magnitude and angle to match power flow result"
    for bus in PSY.get_components(PSY.Bus, sys)
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(inputs)[bus_n]
        initial_guess[bus_ix] = PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
        initial_guess[bus_ix + bus_size] = PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
        CRC.@ignore_derivatives @debug "$(PSY.get_name(bus)) V_r = $(initial_guess[bus_ix]), V_i = $(initial_guess[bus_ix + bus_size])"
    end
    return BUILD_INCOMPLETE
end

function initialize_static_injection!(inputs::SimulationInputs)
    CRC.@ignore_derivatives @debug "Updating Source internal voltage magnitude and angle"
    static_injection_devices = get_static_injectors(inputs)
    parameters = get_parameters(inputs)
    if !isempty(static_injection_devices)
        for s in static_injection_devices
            local_parameters = @view parameters[_get_wrapper_name(s)]
            initialize_static_device!(s, local_parameters)
        end
    end
    return BUILD_INCOMPLETE
end

function initialize_dynamic_injection!(
    initial_guess::AbstractArray,
    inputs::SimulationInputs,
    system::PSY.System,
)
    CRC.@ignore_derivatives @debug "Updating Dynamic Injection Component Initial Guess"
    initial_inner_vars = zeros(get_inner_vars_count(inputs))
    parameters = get_parameters(inputs)
    for dynamic_device in get_dynamic_injectors(inputs)
        static = CRC.@ignore_derivatives PSY.get_component(
            dynamic_device.static_type,
            system,
            PSY.get_name(dynamic_device),
        )
        CRC.@ignore_derivatives @debug "Initializing $(PSY.get_name(dynamic_device)) - $(typeof(dynamic_device.device))"
        _inner_vars = @view initial_inner_vars[get_inner_vars_index(dynamic_device)]
        _parameters = @view parameters[_get_wrapper_name(dynamic_device)]
        _states = @view initial_guess[get_ix_range(dynamic_device)]
        initialize_dynamic_device!(
            dynamic_device,
            static,
            _inner_vars,
            _parameters,
            _states,
        )
    end
    return BUILD_INCOMPLETE
end

function initialize_dynamic_branches!(
    initial_guess::Vector{Float64},
    inputs::SimulationInputs,
)
    parameters = get_parameters(inputs)
    CRC.@ignore_derivatives @debug "Initializing Dynamic Branches"
    for br in get_dynamic_branches(inputs)
        CRC.@ignore_derivatives @debug "$(PSY.get_name(br)) -  $(typeof(br))"
        wrapper_name = _get_wrapper_name(br)
        _parameters = @view parameters[wrapper_name]
        _states = @view initial_guess[get_ix_range(br)]
        initialize_dynamic_device!(br, _parameters, _states)
    end
    return BUILD_INCOMPLETE
end

function check_valid_values(initial_guess::Vector{Float64}, inputs::SimulationInputs)
    invalid_initial_guess = String[]
    for i in get_bus_range(inputs)
        if initial_guess[i] > 1.3 || initial_guess[i] < -1.3
            push!(invalid_initial_guess, "Voltage entry $i")
        end
    end

    for device in get_dynamic_injectors(inputs)
        device_index = get_global_index(device)
        if haskey(device_index, :ω)
            dev_freq = initial_guess[device_index[:ω]]
            if dev_freq > 1.2 || dev_freq < 0.8
                push!(invalid_initial_guess, "$(PSY.get_name(device)) - :ω")
            end
        end
        all(isfinite, initial_guess) && continue
        for state in get_global_index(device)
            if state.second ∈ i
                push!(invalid_initial_guess, "$device - $(state.first)")
            end
        end
    end

    if !isempty(invalid_initial_guess)
        CRC.@ignore_derivatives @error(
            "Invalid initial condition values $invalid_initial_guess"
        )
        return BUILD_FAILED
    end

    return BUILD_IN_PROGRESS
end

#Setter functions for CA
function set_V_ref!(array, value)
    @view(array["refs"])["V_ref"] = value
end
function set_ω_ref!(array, value)
    @view(array["refs"])["ω_ref"] = value
end
function set_θ_ref!(array, value)
    @view(array["refs"])["θ_ref"] = value
end
function set_P_ref!(array, value)
    @view(array["refs"])["P_ref"] = value
end
function set_Q_ref!(array, value)
    @view(array["refs"])["Q_ref"] = value
end

# Default implementation for both models. This implementation is to future proof if there is
# a divergence between the required build methods
#PASS x0 and inputs, but not the full simulation not sim....?
function _initialize_state_space(
    sim::Simulation{T},
    ::Val{POWERFLOW_AND_DEVICES},
) where {T <: SimulationModel}
    inputs = get_simulation_inputs(sim)
    sim.x0 = _get_flat_start(inputs)
    CRC.@ignore_derivatives @info("Pre-Initializing Simulation States")
    @assert sim.status == BUILD_INCOMPLETE
    while sim.status == BUILD_INCOMPLETE
        CRC.@ignore_derivatives @debug "Start state intialization routine"
        sim.status = power_flow_solution!(sim.x0, get_system(sim), inputs)
        sim.status = initialize_static_injection!(inputs)
        sim.status = initialize_dynamic_injection!(sim.x0, inputs, get_system(sim))
        if has_dyn_lines(inputs)
            sim.status = initialize_dynamic_branches!(sim.x0, inputs)
        else
            CRC.@ignore_derivatives @debug "No Dynamic Branches in the system"
        end
        sim.status = check_valid_values(sim.x0, inputs)
    end
end

#GET rid of try catch at that level. 
#GET rid of timing at that level 
function _initialize_state_space(
    sim::Simulation{T},
    ::Val{DEVICES_ONLY},
) where {T <: SimulationModel}
    CRC.@ignore_derivatives @info("Pre-Initializing Simulation States")
    inputs = get_simulation_inputs(sim)
    #GET INNER VARS, PARAMETERS, AND STATES AT THIS LEVEL ... 
    @assert sim.status == BUILD_INCOMPLETE
    while sim.status == BUILD_INCOMPLETE
        CRC.@ignore_derivatives @debug "Start state intialization routine"
        sim.status = initialize_static_injection!(inputs)
        sim.status = initialize_dynamic_injection!(sim.x0, inputs, get_system(sim))
        if has_dyn_lines(inputs)
            sim.status = initialize_dynamic_branches!(sim.x0, inputs)
        else
            CRC.@ignore_derivatives @debug "No Dynamic Branches in the system"
        end
        sim.status = check_valid_values(sim.x0, inputs)
    end
end

function _initialize_state_space(
    sim::Simulation{T},
    ::Val{FLAT_START},
) where {T <: SimulationModel}
    simulation_inputs = get_simulation_inputs(sim)
    sim.x0 = _get_flat_start(simulation_inputs)
end

function _initialize_state_space(
    sim::Simulation{T},
    ::Val{INITIALIZED},
) where {T <: SimulationModel}
    simulation_inputs = get_simulation_inputs(sim)
    @assert sim.status == BUILD_INCOMPLETE
    if length(sim.x0) != get_variable_count(simulation_inputs)
        throw(
            IS.ConflictingInputsError(
                "The size of the provided initial state space does not match the model's state space.",
            ),
        )
    end
    CRC.@ignore_derivatives @warn(
        "Using existing initial conditions value for simulation initialization"
    )
    sim.x0 = deepcopy(sim.x0_init)
    sim.status = SIMULATION_INITIALIZED
    return
end

"""
Returns a Dictionary with the resulting initial conditions of the simulation
"""
function read_initial_conditions(sim::Simulation)
    system = get_system(sim)
    simulation_inputs = get_simulation_inputs(sim)
    bus_size = get_bus_count(simulation_inputs)
    V_R = Dict{Int, Float64}()
    V_I = Dict{Int, Float64}()
    Vm = Dict{Int, Float64}()
    θ = Dict{Int, Float64}()
    for bus in PSY.get_components(PSY.Bus, system)
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(simulation_inputs)[bus_n]
        V_R[bus_n] = get_x0(sim)[bus_ix]
        V_I[bus_n] = get_x0(sim)[bus_ix + bus_size]
        Vm[bus_n] = sqrt(
            get_x0(sim)[bus_ix]^2 +
            get_x0(sim)[bus_ix + bus_size]^2,
        )
        θ[bus_n] = atan(
            get_x0(sim)[bus_ix + bus_size],
            get_x0(sim)[bus_ix],
        )
    end
    results = Dict{String, Any}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in get_dynamic_injectors(simulation_inputs)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = get_global_index(device)
        x0_device = Dict{Symbol, Float64}()
        for s in states
            x0_device[s] = get_x0(sim)[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = get_dynamic_branches(simulation_inputs)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = get_global_index(br)
            x0_br = Dict{Symbol, Float64}()
            for s in states
                x0_br[s] = get_x0(sim)[global_index[s]]
            end
            printed_name = "Line " * name
            results[printed_name] = x0_br
        end
    end
    return results
end

function set_operating_point!(
    x0_init::Vector{Float64},
    inputs::SimulationInputs,
    system::PSY.System,
)
    status = BUILD_INCOMPLETE
    while status == BUILD_INCOMPLETE
        status = power_flow_solution!(x0_init, system, inputs)
        status = initialize_static_injection!(inputs)
        status = initialize_dynamic_injection!(x0_init, inputs, system)
        status = initialize_dynamic_branches!(x0_init, inputs)
        status = SIMULATION_INITIALIZED
    end
    if status != SIMULATION_INITIALIZED
        error("Failed to find operating point")
    end
    return
end
