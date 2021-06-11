function _power_flow_solution!(
    initial_guess::Vector{Float64},
    sys::PSY.System,
    inputs::SimulationInputs,
)
    res = PSY.solve_powerflow!(sys)
    if !res
        @error("PowerFlow failed to solve")
        return BUILD_FAILED
    end
    bus_size = length(PSY.get_bus_numbers(sys))
    @debug "Updating Voltage guess"
    for bus in PSY.get_components(PSY.Bus, sys)
        @debug PSY.get_name(bus)
        @debug "V_r" PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
        @debug "V_i" PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(inputs)[bus_n]
        initial_guess[bus_ix] = PSY.get_magnitude(bus) * cos(PSY.get_angle(bus))
        initial_guess[bus_ix + bus_size] = PSY.get_magnitude(bus) * sin(PSY.get_angle(bus))
    end
    return BUILD_INCOMPLETE
end

function _initialize_static_injection!(inputs::SimulationInputs)
    @debug "Updating Source internal voltages"
    static_injection_devices = get_static_injections_data(inputs)
    if !isempty(static_injection_devices)
        try
            for s in static_injection_devices
                initialize_static_device!(s)
            end
        catch e
            bt = catch_backtrace()
            @error "Static Injection Failed to Initialize" exception = e, bt
            return BUILD_FAILED
        end
    end
    return BUILD_INCOMPLETE
end

function _initialize_dynamic_injection!(
    initial_guess::Vector{Float64},
    inputs::SimulationInputs,
)
    @debug "Updating Dynamic Injection Component Initial Guess"
    injection_start = get_injection_pointer(inputs)
    try
        for d in get_injectors_data(inputs)
            dynamic_device = PSY.get_dynamic_injector(d)
            @debug PSY.get_name(d) typeof(d) typeof(dynamic_device)
            n_states = PSY.get_n_states(dynamic_device)
            ix_range = range(injection_start, length = n_states)
            injection_start = injection_start + n_states
            x0_device = initialize_dynamic_device!(dynamic_device, d)
            @assert length(x0_device) == n_states
            initial_guess[ix_range] = x0_device
        end
    catch e
        bt = catch_backtrace()
        @error "Dynamic Injection Failed to Initialize" exception = e, bt
        return BUILD_FAILED
    end
    return BUILD_INCOMPLETE
end

function _initialize_dynamic_branches!(
    initial_guess::Vector{Float64},
    inputs::SimulationInputs,
)
    @debug "Updating Component Initial Guess"
    branches_start = get_branches_pointer(inputs)
    try
        for br in get_dynamic_branches(inputs)
            @debug PSY.get_name(br) typeof(br)
            n_states = PSY.get_n_states(br)
            ix_range = range(branches_start, length = n_states)
            branches_start = branches_start + n_states
            x0_branch = initialize_dynamic_device!(br)
            @assert length(x0_branch) == n_states
            initial_guess[ix_range] = x0_branch
        end
    catch e
        bt = catch_backtrace()
        @error "Dynamic Branches Failed to Initialize" exception = e, bt
        return BUILD_FAILED
    end
    return BUILD_INCOMPLETE
end

function _check_valid_values(initial_guess::Vector{Float64}, inputs::SimulationInputs)
    if any(!isfinite, initial_guess)
        i = findall(!isfinite, initial_guess)
        invalid_initial_guess = String[]
        for (device, states) in get_global_index(inputs)
            for state in states
                if state.second ∈ i
                    push!(invalid_initial_guess, "$device - $(p.first)")
                end
            end
        end
        @error("Invalid initial guess values $invalid_initial_guess")
        return BUILD_FAILED
    end
    return BUILD_INCOMPLETE
end

struct NLsolveWrapper
    zero::Vector{Float64}
    converged::Bool
    failed::Bool
end

NLsolveWrapper() = NLsolveWrapper(Vector{Float64}(), false, true)
converged(sol::NLsolveWrapper) = sol.converged
failed(sol::NLsolveWrapper) = sol.failed

function _nlsolve_call(initial_guess, inputs::SimulationInputs, tolerance::Float64)
    dx0 = zeros(length(initial_guess)) #Define a vector of zeros for the derivative
    inif! = (out, x) -> system_implicit!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        inputs,    #Parameters
        -99.0,    #time val not relevant
    )
    return try
        sys_solve = NLsolve.nlsolve(
            inif!,
            initial_guess,
            xtol = tolerance,
            ftol = tolerance,
            method = :trust_region,
            #autodiff = :forward,
        ) #Solve using initial guess x0
        NLsolveWrapper(sys_solve.zero, NLsolve.converged(sys_solve), false)
    catch e
        bt = catch_backtrace()
        @error "NLsolve failed to solve" exception = e, bt
        NLsolveWrapper()
    end
end

function _convergence_check(sys_solve::NLsolveWrapper, tol::Float64)
    if converged(sys_solve)
        @info("Initialization succeeded with a tolerance of $(tol). Saving solution")
    else
        @warn(
            "Initialization convergence failed, initial conditions do not meet conditions for an stable equilibrium.\nTrying to solve again reducing numeric tolerance"
        )
    end
    return converged(sys_solve)
end

function _refine_initial_condition!(
    initial_guess::Vector{Float64},
    inputs::SimulationInputs,
)
    @debug "Start NLSolve System Run"
    @debug initial_guess
    converged = false
    for tol in [STRICT_NL_SOLVE_TOLERANCE, RELAXED_NL_SOLVE_TOLERANCE]
        sys_solve = _nlsolve_call(initial_guess, inputs, tol)
        failed(sys_solve) && return BUILD_FAILED
        converged = _convergence_check(sys_solve, tol)
        @debug "Write result to initial guess vector under condition converged = $(converged)"
        initial_guess .= sys_solve.zero
        if converged
            break
        end
        @warn(
            "Initialization never converged to desired tolerances. Initial conditions do not meet conditions for an stable equilibrium. Simulation migth diverge"
        )
    end
    return converged ? SIMULATION_INITIALIZED : _check_valid_values(initial_guess, inputs)
end

# Default implementation for both models. This implementation is to future proof if there is
# a divergence between the required build methods
function _calculate_initial_conditions!(sim::Simulation)
    inputs = get_simulation_inputs(sim)
    simulation_pre_step!(inputs, Real)
    @debug "Start state intialization routine"
    while sim.status == BUILD_INCOMPLETE
        sim.status = _power_flow_solution!(sim.x0_init, get_system(sim), inputs)
        sim.status = _initialize_static_injection!(inputs)
        sim.status = _initialize_dynamic_injection!(sim.x0_init, inputs)
        if has_dyn_lines(inputs)
            sim.status = _initialize_dynamic_branches!(sim.x0_init, inputs)
        else
            @debug "No Dynamic Branches in the system"
        end
        sim.status = _check_valid_values(sim.x0_init, inputs)
        sim.status = _refine_initial_condition!(sim.x0_init, inputs)
    end
    return sim.status != BUILD_FAILED
end

function calculate_initial_conditions!(sim::Simulation{ImplicitModel})
    return _calculate_initial_conditions!(sim)
end

function calculate_initial_conditions!(sim::Simulation{MassMatrixModel})
    return _calculate_initial_conditions!(sim)
end

"""
Returns a Dictionary with the resulting initial conditions of the simulation
"""
function get_initial_conditions(sim::Simulation)
    system = get_system(sim)
    bus_size = get_bus_count(sim.inputs)
    V_R = Dict{Int, Float64}()
    V_I = Dict{Int, Float64}()
    Vm = Dict{Int, Float64}()
    θ = Dict{Int, Float64}()
    for bus in PSY.get_components(PSY.Bus, system)
        bus_n = PSY.get_number(bus)
        bus_ix = get_lookup(sim.inputs)[bus_n]
        V_R[bus_n] = sim.x0_init[bus_ix]
        V_I[bus_n] = sim.x0_init[bus_ix + bus_size]
        Vm[bus_n] = sqrt(sim.x0_init[bus_ix]^2 + sim.x0_init[bus_ix + bus_size]^2)
        θ[bus_n] = angle(sim.x0_init[bus_ix] + sim.x0_init[bus_ix + bus_size] * 1im)
    end
    results = Dict{String, Any}("V_R" => V_R, "V_I" => V_I, "Vm" => Vm, "θ" => θ)
    for device in PSY.get_components(PSY.DynamicInjection, system)
        states = PSY.get_states(device)
        name = PSY.get_name(device)
        global_index = get_global_index(sim.inputs)[name]
        x0_device = Dict{Symbol, Float64}()
        for s in states
            x0_device[s] = sim.x0_init[global_index[s]]
        end
        results[name] = x0_device
    end
    dyn_branches = PSY.get_components(PSY.DynamicBranch, system)
    if !isempty(dyn_branches)
        for br in dyn_branches
            states = PSY.get_states(br)
            name = PSY.get_name(br)
            global_index = get_global_index(sim.inputs)[name]
            x0_br = Dict{Symbol, Float64}()
            for s in states
                x0_br[s] = sim.x0_init[global_index[s]]
            end
            printed_name = "Line " * name
            results[printed_name] = x0_br
        end
    end
    return results
end
