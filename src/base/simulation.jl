mutable struct Simulation{T <: SimulationModel}
    status::BUILD_STATUS
    problem::Union{Nothing, SciMLBase.DEProblem}
    sys::PSY.System
    perturbations::Vector{<:Perturbation}
    x0_init::Vector{Float64}
    initialized::Bool
    tstops::Vector{Float64}
    callbacks::DiffEqBase.CallbackSet
    solution::Union{Nothing, SciMLBase.AbstractODESolution}
    simulation_folder::String
    simulation_inputs::SimulationInputs
    console_level::Base.CoreLogging.LogLevel
    file_level::Base.CoreLogging.LogLevel
    multimachine::Bool
end

get_system(sim::Simulation) = sim.sys

function Simulation(
    ::Type{T},
    sys::PSY.System;
    simulation_inputs,
    perturbations = Vector{Perturbation}(),
    simulation_folder::String = "",
    console_level = Logging.Warn,
    file_level = Logging.Debug,
) where {T <: SimulationModel}
    return Simulation{T}(
        BUILD_INCOMPLETE,
        nothing,
        sys,
        perturbations,
        Vector{Float64}(),
        false,
        Vector{Float64}(),
        DiffEqBase.CallbackSet(),
        nothing,
        simulation_folder,
        simulation_inputs,
        console_level,
        file_level,
        false,
    )
end

function Simulation!(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    tspan::NTuple{2, Float64},
    perturbation::Perturbation;
    kwargs...,
) where {T <: SimulationModel}
    return Simulation!(T, system, simulation_folder, tspan, [perturbation]; kwargs...)
end

function Simulation(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    tspan::NTuple{2, Float64},
    perturbation::Perturbation;
    kwargs...,
) where {T <: SimulationModel}
    return Simulation(T, system, simulation_folder, tspan, [perturbation]; kwargs...)
end

"""
Builds the simulation object and conducts the indexing process. The initial conditions are stored in the system.

# Accepted Key Words
- `initialize_simulation::Bool : Runs the initialization routine. If false, simulation runs based on the operation point stored in System`
- `system_to_file::Bool`: Serializes the initialized system
- `console_level::Logging`: Sets the level of logging output to the console. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debug
- `file_level::Logging`: Sets the level of logging output to file. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debug
"""
function Simulation!(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
) where {T <: SimulationModel}
    check_folder(simulation_folder)
    # Instantiates the Simulation object
    sim = Simulation(
        T,
        system;
        simulation_inputs = SimulationInputs(sys = system, tspan = tspan),
        simulation_folder = simulation_folder,
        perturbations = perturbations,
        console_level = get(kwargs, :console_level, Logging.Warn),
        file_level = get(kwargs, :file_level, Logging.Debug),
    )
    build!(sim; kwargs...)
    if get(kwargs, :system_to_file, false)
        PSY.to_json(system, joinpath(simulation_folder, "initialized_system.json"))
    end
    return sim
end

"""
Initializes the simulations and builds the indexing. The input system is not modified during the initialization

# Accepted Key Words
- `initialize_simulation::Bool : Runs the initialization routine. If false, simulation runs based on the operation point stored in System`
- `system_to_file::Bool`: Serializes the original input system
- `console_level::Logging`: Sets the level of logging output to the console. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debugg
- `file_level::Logging`: Sets the level of logging output to file. Can be set to Logging.Error, Logging.Warn, Logging.Info or Logging.Debug
"""
function Simulation(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
) where {T <: SimulationModel}
    check_folder(simulation_folder)
    simulation_system = deepcopy(system)
    # Instantiates the Simulation object
    sim = Simulation(
        T,
        simulation_system;
        simulation_inputs = SimulationInputs(sys = simulation_system, tspan = tspan),
        simulation_folder = simulation_folder,
        perturbations = perturbations,
        console_level = get(kwargs, :console_level, Logging.Warn),
        file_level = get(kwargs, :file_level, Logging.Debug),
    )
    build!(sim; kwargs...)
    if get(kwargs, :system_to_file, false)
        PSY.to_json(system, joinpath(simulation_folder, "input_system.json"))
    end
    return sim
end

function reset!(sim::Simulation{T}) where {T <: SimulationModel}
    @info "Rebuilding the simulation after reset"
    sim.simulation_inputs =
        SimulationInputs(sys = get_system(sim), tspan = sim.simulation_inputs.tspan)
    build!(sim; file_mode = "a")
    @info "Simulation reset to status $(sim.status)"
    return
end

function configure_logging(sim::Simulation, file_mode)
    return IS.configure_logging(
        console = true,
        console_stream = stderr,
        console_level = sim.console_level,
        file = true,
        filename = joinpath(sim.simulation_folder, SIMULATION_LOG_FILENAME),
        file_level = sim.file_level,
        file_mode = file_mode,
        tracker = nothing,
        set_global = false,
    )
end

function build!(sim; file_mode = "w", kwargs...)
    logger = configure_logging(sim, file_mode)
    try
        Logging.with_logger(logger) do
            _build!(sim; kwargs...)
        end
    finally
        close(logger)
    end
    return
end

function _build!(sim::Simulation{T}; kwargs...) where {T <: SimulationModel}
    simulation_system = get_system(sim)
    sim.status = BUILD_INCOMPLETE
    PSY.set_units_base_system!(simulation_system, "DEVICE_BASE")
    check_kwargs(kwargs, SIMULATION_ACCEPTED_KWARGS, "Simulation")
    simulation_inputs = build!(sim.simulation_inputs, T, simulation_system)
    @debug "Simulation Inputs Created"
    var_count = get_variable_count(simulation_inputs)

    flat_start = zeros(var_count)
    bus_count = length(PSY.get_components(PSY.Bus, simulation_system))
    flat_start[1:bus_count] .= 1.0
    sim.x0_init = get(kwargs, :initial_guess, flat_start)

    initialize_simulation = get(kwargs, :initialize_simulation, true)
    if initialize_simulation
        @info("Initializing Simulation States")
        add_aux_arrays!(simulation_inputs, Float64)
        sim.initialized = calculate_initial_conditions!(sim, simulation_inputs)
        if sim.status == BUILD_FAILED
            @error("Simulation Build Failed. Simulations status = $(sim.status)")
            return
        end
    end

    dx0 = zeros(var_count)
    _build_perturbations!(sim)
    add_aux_arrays!(simulation_inputs, Float64)
    sim.problem = SciMLBase.DAEProblem(
        system!,
        dx0,
        sim.x0_init,
        get_tspan(sim.simulation_inputs),
        simulation_inputs,
        differential_vars = get_DAE_vector(simulation_inputs);
        kwargs...,
    )
    sim.multimachine = (get_global_vars(simulation_inputs)[:ω_sys_index] != 0)
    sim.status = BUILT
    @info "Completed Build Successfully. Simulations status = $(sim.status)"
    return
end

function _build_perturbations!(sim::Simulation)
    @info "Attaching Perturbations"
    if isempty(sim.perturbations)
        @debug "The simulation has no perturbations"
        return DiffEqBase.CallbackSet(), [0.0]
    end
    system = get_system(sim)
    perturbations = sim.perturbations
    perturbations_count = length(perturbations)
    callback_vector = Vector{DiffEqBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Vector{Float64}(undef, perturbations_count)
    for (ix, pert) in enumerate(perturbations)
        @debug pert
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(system, pert)
        callback_vector[ix] = DiffEqBase.DiscreteCallback(condition, affect)
        tstops[ix] = pert.time
    end
    sim.tstops = tstops
    sim.callbacks = DiffEqBase.CallbackSet((), tuple(callback_vector...))
    return
end

function _simulation_pre_step(sim::Simulation, reset_simulation::Bool)
    if sim.status != BUILT && !reset_simulation
        error(
            "The Simulation status is $(sim.status). Use keyword argument reset_simulation = true",
        )
    end

    if reset_simulation
        reset!(sim)
    end
    return
end

function execute!(sim::Simulation, solver; kwargs...)
    @debug "status before execute" sim.status
    reset_simulation = get(kwargs, :reset_simulation, false)
    reset_simulation = sim.status == CONVERTED_FOR_SMALL_SIGNAL || reset_simulation
    _simulation_pre_step(sim, reset_simulation)
    sim.status = SIMULATION_STARTED
    sim.solution = SciMLBase.solve(
        sim.problem,
        solver;
        callback = sim.callbacks,
        tstops = sim.tstops,
        kwargs...,
    )
    if sim.solution.retcode == :Success
        sim.status = SIMULATION_FINALIZED
    else
        sim.status = SIMULATION_FAILED
    end
    return
end
