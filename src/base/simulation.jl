mutable struct Simulation{T <: SimulationModel}
    status::BUILD_STATUS
    problem::Union{Nothing, SciMLBase.DEProblem}
    tspan::NTuple{2, Float64}
    sys::PSY.System
    perturbations::Vector{<:Perturbation}
    x0_init::Vector{Float64}
    initialized::Bool
    tstops::Vector{Float64}
    callbacks::DiffEqBase.CallbackSet
    simulation_folder::String
    inputs::Union{Nothing, SimulationInputs}
    results::Union{Nothing, SimulationResults}
    console_level::Base.CoreLogging.LogLevel
    file_level::Base.CoreLogging.LogLevel
    multimachine::Bool
end

get_system(sim::Simulation) = sim.sys
get_simulation_inputs(sim::Simulation) = sim.inputs
get_initial_conditions(sim::Simulation) = sim.x0_init
get_tspan(sim::Simulation) = sim.tspan

function Simulation(
    ::Type{T},
    sys::PSY.System;
    tspan,
    initial_conditions,
    initialize_simulation,
    perturbations,
    simulation_folder,
    console_level,
    file_level,
) where {T <: SimulationModel}
    PSY.set_units_base_system!(sys, "DEVICE_BASE")

    return Simulation{T}(
        BUILD_INCOMPLETE,
        nothing,
        tspan,
        sys,
        perturbations,
        initial_conditions,
        !initialize_simulation,
        Vector{Float64}(),
        DiffEqBase.CallbackSet(),
        simulation_folder,
        nothing,
        nothing,
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
- `initial_conditions::Vector{Float64} : Allows the user to pass a vector with the initial condition values desired in the simulation. If initialize_simulation = true, these values are used as a first guess and overwritten.
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
        tspan = tspan,
        initial_conditions = get(kwargs, :initial_conditions, Vector{Float64}()),
        initialize_simulation = get(kwargs, :initialize_simulation, true),
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
- `initial_conditions::Vector{Float64} : Initial conditions provided by the user. Must be used with initialize_simulation = false to not be override.`
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
        tspan = tspan,
        initial_conditions = get(kwargs, :initial_conditions, Vector{Float64}()),
        initialize_simulation = get(kwargs, :initialize_simulation, true),
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
    sim.inputs = SimulationInputs(T(), get_system(sim), sim.inputs.tspan)
    build!(sim)
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

function _build_inputs!(sim::Simulation{T}) where {T <: SimulationModel}
    simulation_system = get_system(sim)
    sim.inputs = SimulationInputs(T, simulation_system)
    @debug "Simulation Inputs Created"
    return
end

function _get_flat_start(inputs::SimulationInputs)
    bus_count = get_bus_count(inputs)
    var_count = get_variable_count(inputs)
    initial_conditions = zeros(var_count)
    initial_conditions[1:bus_count] .= 1.0
    return initial_conditions
end

function _initialize_state_space(sim::Simulation{T}) where {T <: SimulationModel}
    simulation_inputs = get_simulation_inputs(sim)
    x0_init = sim.x0_init
    if isempty(x0_init) && sim.initialized
        @warn "Initial Conditions set to flat start"
        sim.x0_init = _get_flat_start(simulation_inputs)
    elseif isempty(x0_init) && !sim.initialized
        sim.x0_init = _get_flat_start(simulation_inputs)
    elseif !isempty(x0_init) && sim.initialized
        if length(sim.x0_init) != get_variable_count(simulation_inputs)
            throw(
                IS.ConflictingInputsError(
                    "The size of the provided initial state space does not match the model's state space.",
                ),
            )
        end
    elseif !isempty(x0_init) && !sim.initialized
        @warn "initial_conditions were provided with initialize_simulation. User's initial_conditions will be overwritten."
        if length(sim.x0_init) != get_variable_count(simulation_inputs)
            sim.x0_init = _get_flat_start(simulation_inputs)
        end
    else
        @assert false
    end
    return
end

function _pre_initialize_simulation!(sim::Simulation)
    _initialize_state_space(sim)
    if sim.initialized != true
        @info("Pre-Initializing Simulation States")
        sim.initialized = precalculate_initial_conditions!(sim)
        if !sim.initialized
            error(
                "The simulation failed to find an adequate initial guess for the initialization. Check the intialization routine.",
            )
        end
    else
        @warn(
            "No Pre-initialization conducted. If this is unexpected, check the initialization keywords"
        )
        sim.status = SIMULATION_INITIALIZED
    end
    return
end

function _get_jacobian(sim::Simulation{T}) where {T <: SimulationModel}
    inputs = get_simulation_inputs(sim)
    x0_init = get_initial_conditions(sim)
    return JacobianFunctionWrapper(T(inputs, x0_init, JacobianCache), x0_init)
end

function _build_perturbations!(sim::Simulation)
    @info "Attaching Perturbations"
    if isempty(sim.perturbations)
        @debug "The simulation has no perturbations"
        return DiffEqBase.CallbackSet(), [0.0]
    end
    inputs = get_simulation_inputs(sim)
    perturbations = sim.perturbations
    perturbations_count = length(perturbations)
    callback_vector = Vector{DiffEqBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Vector{Float64}(undef, perturbations_count)
    for (ix, pert) in enumerate(perturbations)
        @debug pert
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(inputs, get_system(sim), pert)
        callback_vector[ix] = DiffEqBase.DiscreteCallback(condition, affect)
        tstops[ix] = pert.time
    end
    sim.tstops = tstops
    sim.callbacks = DiffEqBase.CallbackSet((), tuple(callback_vector...))
    return
end

function _get_diffeq_problem(
    sim::Simulation,
    model::SystemModel{ResidualModel},
    jacobian::JacobianFunctionWrapper,
)
    x0 = get_initial_conditions(sim)
    dx0 = zeros(length(x0))
    simulation_inputs = get_simulation_inputs(sim)
    sim.problem = SciMLBase.DAEProblem(
        SciMLBase.DAEFunction{true}(
            model;
            # Currently commented for Sundials compatibility
            #jac = jacobian,
            tgrad = (dT, u, p, t) -> dT .= false,
            #jac_prototype = jacobian.Jv,
        ),
        dx0,
        x0,
        get_tspan(sim),
        simulation_inputs,
        differential_vars = get_DAE_vector(simulation_inputs),
    )
    sim.status = BUILT
    return
end

function _get_diffeq_problem(
    sim::Simulation,
    model::SystemModel{MassMatrixModel},
    jacobian::JacobianFunctionWrapper,
)
    simulation_inputs = get_simulation_inputs(sim)
    sim.problem = SciMLBase.ODEProblem(
        SciMLBase.ODEFunction{true}(
            model,
            mass_matrix = get_mass_matrix(simulation_inputs),
            jac = jacobian,
            jac_prototype = jacobian.Jv,
            # Necessary to avoid unnecessary calculations in Rosenbrock methods
            tgrad = (dT, u, p, t) -> dT .= false,
        ),
        sim.x0_init,
        get_tspan(sim),
        simulation_inputs,
    )
    sim.status = BUILT
    return
end

function _build!(sim::Simulation{T}; kwargs...) where {T <: SimulationModel}
    check_kwargs(kwargs, SIMULATION_ACCEPTED_KWARGS, "Simulation")
    _build_inputs!(sim)
    _pre_initialize_simulation!(sim)
    if sim.status != BUILD_FAILED
        simulation_inputs = get_simulation_inputs(sim)
        try
            jacobian = _get_jacobian(sim)
            model = T(simulation_inputs, get_initial_conditions(sim), SimCache)
            refine_initial_condition!(sim, model, jacobian)
            _build_perturbations!(sim)
            _get_diffeq_problem(sim, model, jacobian)
            @info "Simulations status = $(sim.status)"
        catch e
            bt = catch_backtrace()
            @error "$T failed to build" exception = e, bt
            sim.status = BUILD_FAILED
        end
    else
        @error "The simulation couldn't be initialized correctly. Simulations status = $(sim.status)"
    end
    return
end

function build!(sim; kwargs...)
    logger = configure_logging(sim, "w")
    try
        Logging.with_logger(logger) do
            _build!(sim; kwargs...)
        end
    finally
        close(logger)
    end
    return
end

function simulation_pre_step!(sim::Simulation, reset_sim::Bool, ::Type{T}) where {T <: Real}
    reset_sim = sim.status == CONVERTED_FOR_SMALL_SIGNAL || reset_sim
    if sim.status == BUILD_FAILED
        error(
            "The Simulation status is $(sim.status). Can not continue, correct your inputs and build the simulation again.",
        )
    elseif sim.status != BUILT && !reset_sim
        error(
            "The Simulation status is $(sim.status). Use keyword argument reset_simulation = true",
        )
    end

    reset_sim && reset!(sim)
    return
end

function execute!(sim::Simulation, solver; kwargs...)
    @debug "status before execute" sim.status
    simulation_pre_step!(sim, get(kwargs, :reset_simulation, false), Float64)
    sim.status = SIMULATION_STARTED
    time_log = Dict{Symbol, Any}()
    solution,
    time_log[:timed_solve_time],
    time_log[:solve_bytes_alloc],
    time_log[:sec_in_gc] = @timed SciMLBase.solve(
        sim.problem,
        solver;
        callback = sim.callbacks,
        tstops = sim.tstops,
        kwargs...,
    )
    if solution.retcode == :Success
        sim.status = SIMULATION_FINALIZED
        sim.results = SimulationResults(
            get_simulation_inputs(sim),
            get_system(sim),
            time_log,
            solution,
        )
    else
        @error("The simulation failed with return code $(solution.retcode)")
        sim.status = SIMULATION_FAILED
    end
    return sim.status
end

function read_results(sim::Simulation)
    return sim.results
end
