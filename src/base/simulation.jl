mutable struct Simulation{T <: SimulationModel}
    status::BUILD_STATUS
    sys::PSY.System
    perturbations::Vector{<:Perturbation}
    simulation_folder::String
    inputs::Union{Nothing, SimulationInputs}
    results::Union{Nothing, SimulationResults}
    internal::Union{Nothing, SimulationInternal}
    console_level::Base.CoreLogging.LogLevel
    file_level::Base.CoreLogging.LogLevel
    multimachine::Bool
    user_provided_initial_conditions::Vector{Float64}
    initialized::Bool
end

get_system(sim::Simulation) = sim.sys
get_simulation_inputs(sim::Simulation) = sim.inputs
get_initial_conditions(sim::Simulation) = sim.internal.x0_init
get_tspan(sim::Simulation) = sim.internal.tspan

function Simulation(
    ::Type{T},
    sys::PSY.System;
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
        sys,
        perturbations,
        simulation_folder,
        nothing,
        nothing,
        nothing,
        console_level,
        file_level,
        false,
        initial_conditions,
        !initialize_simulation,
    )
end

function Simulation!(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    perturbation::Perturbation;
    kwargs...,
) where {T <: SimulationModel}
    return Simulation!(T, system, simulation_folder, [perturbation]; kwargs...)
end

function Simulation(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    perturbation::Perturbation;
    kwargs...,
) where {T <: SimulationModel}
    return Simulation(T, system, simulation_folder, [perturbation]; kwargs...)
end

"""
    function Simulation!
        ::SimulationModel
        system::PowerSystems.System
        simulation_folder::String
        perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
        kwargs...,
    end

Builds the simulation object and conducts the indexing process. The original system is not modified and a copy its created and stored in the Simulation.

# Arguments:
- `::SimulationModel` : Type of Simulation Model. `ResidualModel` or `MassMatrixModel`. See [Models Section](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/models/) for more details
- `system::PowerSystems.System` : System data
- `simulation_folder::String` : Folder directory
- `tspan::NTuple{2, Float64}` : Time span for simulation
- `perturbations::Vector{<:Perturbation}` : Vector of Perturbations for the Simulation. Default: No Perturbations
- `initialize_simulation::Bool` : Runs the initialization routine. If false, simulation runs based on the operation point stored in System
- `initial_conditions::Vector{Float64}` : Allows the user to pass a vector with the initial condition values desired in the simulation. If initialize_simulation = true, these values are used as a first guess and overwritten.
- `frequency_reference` : Default `ReferenceBus`. Determines which frequency model is used for the network. Currently there are two options available:
    - `ConstantFrequency` assumes that the network frequency is 1.0 per unit at all times.
    - `ReferenceBus` will use the frequency state of a Dynamic Generator (rotor speed) or Dynamic Inverter (virtual speed) connected to the Reference Bus (defined in the Power Flow data) as the network frequency. If multiple devices are connected to such bus, the device with larger base power will be used as a reference. If a Voltage Source is connected to the Reference Bus, then a `ConstantFrequency` model will be used.
- `system_to_file::Bool` : Default `false`. Serializes the initialized system
- `console_level::Logging` : Default `Logging.Warn`. Sets the level of logging output to the console. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`
- `file_level::Logging` : Default `Logging.Debug`. Sets the level of logging output to file. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`
- `disable_timer_output::Bool` : Default `false`. Allows the user to display timer information about the construction and initilization of the Simulation.
"""
function Simulation!(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
) where {T <: SimulationModel}
    check_folder(simulation_folder)
    # Instantiates the Simulation object
    sim = Simulation(
        T,
        system;
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
    function Simulation
        ::SimulationModel
        system::PowerSystems.System
        simulation_folder::String
        perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
        kwargs...,
    end

Builds the simulation object and conducts the indexing process. The initial conditions are stored in the system.

# Arguments:
- `::SimulationModel` : Type of Simulation Model. `ResidualModel` or `MassMatrixModel`. See [Models Section](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/models/) for more details
- `system::PowerSystems.System` : System data
- `simulation_folder::String` : Folder directory
- `tspan::NTuple{2, Float64}` : Time span for simulation
- `perturbations::Vector{<:Perturbation}` : Vector of Perturbations for the Simulation. Default: No Perturbations
- `initialize_simulation::Bool` : Runs the initialization routine. If false, simulation runs based on the operation point stored in System
- `initial_conditions::Vector{Float64}` : Allows the user to pass a vector with the initial condition values desired in the simulation. If initialize_simulation = true, these values are used as a first guess and overwritten.
- `frequency_reference` : Default `ReferenceBus`. Determines which frequency model is used for the network. Currently there are two options available:
    - `ConstantFrequency` assumes that the network frequency is 1.0 per unit at all times.
    - `ReferenceBus` will use the frequency state of a Dynamic Generator (rotor speed) or Dynamic Inverter (virtual speed) connected to the Reference Bus (defined in the Power Flow data) as the network frequency. If multiple devices are connected to such bus, the device with larger base power will be used as a reference. If a Voltage Source is connected to the Reference Bus, then a `ConstantFrequency` model will be used.
- `system_to_file::Bool` : Default `false`. Serializes the initialized system
- `console_level::Logging` : Default `Logging.Warn`. Sets the level of logging output to the console. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`
- `file_level::Logging` : Default `Logging.Debug`. Sets the level of logging output to file. Can be set to `Logging.Error`, `Logging.Warn`, `Logging.Info` or `Logging.Debug`
- `disable_timer_output::Bool` : Default `false`. Allows the user to display timer information about the construction and initilization of the Simulation.
"""
function Simulation(
    ::Type{T},
    system::PSY.System,
    simulation_folder::String,
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
) where {T <: SimulationModel}
    check_folder(simulation_folder)
    simulation_system = deepcopy(system)
    # Instantiates the Simulation object
    sim = Simulation(
        T,
        simulation_system;
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

function configure_logging(sim::Simulation, file_mode; kwargs...)
    return IS.configure_logging(
        console = true,
        console_stream = stderr,
        console_level = get(kwargs, :console_level, sim.console_level),
        file = true,
        filename = joinpath(sim.simulation_folder, SIMULATION_LOG_FILENAME),
        file_level = get(kwargs, :file_level, sim.file_level),
        file_mode = file_mode,
        tracker = nothing,
        set_global = false,
        progress = _prog_meter_enabled(),
    )
end

function _build_inputs!(
    sim::Simulation{T},
    frequency_reference,
) where {T <: SimulationModel}
    simulation_system = get_system(sim)
    sim.inputs = SimulationInputs(T, simulation_system, frequency_reference)
    @debug "Simulation Inputs Created"
    return
end

function _initialize_state_space(sim::Simulation{T}) where {T <: SimulationModel}
    simulation_inputs = get_simulation_inputs(sim)
    x0_init = deepcopy(sim.user_provided_initial_conditions)
    if isempty(x0_init) && sim.initialized
        @warn "Initial Conditions set to flat start"
        x0_init = get_flat_start(simulation_inputs)
    elseif isempty(x0_init) && !sim.initialized
        x0_init = get_flat_start(simulation_inputs)
    elseif !isempty(x0_init) && sim.initialized
        if length(x0_init) != get_variable_count(simulation_inputs)
            throw(
                IS.ConflictingInputsError(
                    "The size of the provided initial state space does not match the model's state space.",
                ),
            )
        end
    elseif !isempty(x0_init) && !sim.initialized
        @warn "initial_conditions were provided with initialize_simulation. User's initial_conditions will be overwritten."
        if length(x0_init) != get_variable_count(simulation_inputs)
            x0_init = get_flat_start(simulation_inputs)
        end
    else
        @assert false
    end
    return x0_init
end

function _pre_initialize_simulation!(sim::Simulation)
    x0_init = _initialize_state_space(sim)
    if sim.initialized != true
        @info("Pre-Initializing Simulation States")
        sim.initialized = precalculate_initial_conditions!(x0_init, sim)
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
    return x0_init
end

function _get_jacobian(::Simulation{T}, inputs, x0_init) where {T <: SimulationModel}
    return get_jacobian(T, inputs, x0_init)
end

function _build_perturbations!(sim::Simulation)
    @info "Attaching Perturbations"
    if isempty(sim.perturbations)
        @debug "The simulation has no perturbations"
        return SciMLBase.CallbackSet(), [0.0]
    end
    inputs = get_simulation_inputs(sim)
    perturbations = sim.perturbations
    perturbations_count = length(perturbations)
    callback_vector = Vector{SciMLBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Vector{Float64}(undef, perturbations_count)
    for (ix, pert) in enumerate(perturbations)
        @debug pert
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(inputs, get_system(sim), pert)
        callback_vector[ix] = SciMLBase.DiscreteCallback(condition, affect)
        tstops[ix] = pert.time
    end
    sim.internal.tstops = tstops
    sim.internal.callbacks = SciMLBase.CallbackSet((), tuple(callback_vector...))
    return
end

function _get_model_function(
    ::Simulation,
    model::SystemModel{ResidualModel},
    jacobian::JacobianFunctionWrapper,
)
    return SciMLBase.DAEFunction{true}(
        model;
        # Currently commented for Sundials.jl compatibility
        #jac = jacobian,
        tgrad = (dT, u, p, t) -> dT .= false,
        #jac_prototype = jacobian.Jv,
    )
end

function _get_model_function(
    sim::Simulation,
    model::SystemModel{MassMatrixModel},
    jacobian::JacobianFunctionWrapper,
)
    simulation_inputs = get_simulation_inputs(sim)
    return SciMLBase.ODEFunction{true}(
        model,
        mass_matrix = get_mass_matrix(simulation_inputs),
        jac = jacobian,
        jac_prototype = jacobian.Jv,
        # Necessary to avoid unnecessary calculations in Rosenbrock methods
        tgrad = (dT, u, p, t) -> dT .= false,
    )
end

function _build!(sim::Simulation{T}; kwargs...) where {T <: SimulationModel}
    check_kwargs(kwargs, SIMULATION_ACCEPTED_KWARGS, "Simulation")
    # Branches are a super set of Lines. Passing both kwargs will
    # be redundant.
    TimerOutputs.reset_timer!(BUILD_TIMER)
    if get(kwargs, :disable_timer_outputs, false)
        TimerOutputs.disable_timer!(BUILD_PROBLEMS_TIMER)
    end

    TimerOutputs.@timeit BUILD_TIMER "Build Simulation" begin
        if get(kwargs, :all_branches_dynamic, false)
            TimerOutputs.@timeit BUILD_TIMER "AC Branch Transform to Dynamic" begin
                sys = get_system(sim)
                transform_branches_to_dynamic(sys, PSY.ACBranch)
            end
        elseif get(kwargs, :all_lines_dynamic, false)
            TimerOutputs.@timeit BUILD_TIMER "Line Transform to Dynamic" begin
                sys = get_system(sim)
                transform_branches_to_dynamic(sys, PSY.Line)
            end
        end
        TimerOutputs.@timeit BUILD_TIMER "Build Simulation Inputs" begin
            f_ref = get(kwargs, :frequency_reference, ReferenceBus)
            _build_inputs!(sim, f_ref)
            # TODO: Update and store f_ref somewhere.
            sim.multimachine =
                get_global_vars_update_pointers(sim.inputs)[GLOBAL_VAR_SYS_FREQ_INDEX] != 0
        end
        TimerOutputs.@timeit BUILD_TIMER "Pre-initialization" begin
            x0_init = _pre_initialize_simulation!(sim)
        end
        if sim.status != BUILD_FAILED
            simulation_inputs = get_simulation_inputs(sim)
            try
                TimerOutputs.@timeit BUILD_TIMER "Calculate Jacobian" begin
                    jacobian = _get_jacobian(sim, simulation_inputs, x0_init)
                end
                TimerOutputs.@timeit BUILD_TIMER "Make Model Function" begin
                    model = T(simulation_inputs, x0_init, SimCache)
                end
                TimerOutputs.@timeit BUILD_TIMER "Initial Condition NLsolve refinement" begin
                    refine_initial_condition!(x0_init, sim, model, jacobian)
                end
                TimerOutputs.@timeit BUILD_TIMER "Make DiffEq Problem" begin
                    model_function = _get_model_function(sim, model, jacobian)
                end

                sim.internal = SimulationInternal(x0_init, jacobian, model, model_function)
                TimerOutputs.@timeit BUILD_TIMER "Build Perturbations" begin
                    _build_perturbations!(sim)
                end
                sim.status = BUILT
                @info "Simulations status = $(sim.status)"
            catch e
                bt = catch_backtrace()
                @error "$T failed to build" exception = e, bt
                sim.status = BUILD_FAILED
            end
        else
            @error "The simulation couldn't be initialized correctly. Simulations status = $(sim.status)"
        end
    end
    return
end

function build!(sim; kwargs...)
    logger = configure_logging(sim, "w")
    try
        Logging.with_logger(logger) do
            _build!(sim; kwargs...)
        end
    catch e
        @error "Build failed" exception = (e, catch_backtrace())
        return sim.status
    finally
        string_buffer = IOBuffer()
        TimerOutputs.print_timer(
            string_buffer,
            BUILD_TIMER,
            sortby = :firstexec,
            compact = true,
        )
        @info "\n$(String(take!(string_buffer)))\n"
        close(logger)
    end
    return
end

function simulation_pre_step!(sim::Simulation, reset_sim::Bool)
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

function _prog_meter_enabled()
    return isa(stderr, Base.TTY) &&
           (get(ENV, "CI", nothing) != "true") &&
           (get(ENV, "RUNNING_PSID_TESTS", nothing) != "true")
end

function _get_diff_eq_problem(
    internal::SimulationInternal{T, SystemModel{ResidualModel, W}, V},
    simulation_inputs::SimulationInputs,
    tspan,
) where {T, V, W}
    x0 = get_initial_conditions(internal)
    dx0 = similar(x0)
    return SciMLBase.DAEProblem(
        internal.sciml_function,
        dx0,
        x0,
        tspan,
        simulation_inputs,
        differential_vars = get_DAE_vector(simulation_inputs),
    )
end

function _get_diff_eq_problem(
    internal::SimulationInternal{T, SystemModel{MassMatrixModel, W}, V},
    simulation_inputs::SimulationInputs,
    tspan,
) where {T, V, W}
    x0 = get_initial_conditions(internal)
    return SciMLBase.ODEProblem(internal.sciml_function, x0, tspan, simulation_inputs)
end

function _execute!(sim::Simulation, solver, tspan; kwargs...)
    @debug "status before execute" sim.status
    simulation_pre_step!(sim, get(kwargs, :reset_simulation, false))
    sim.status = SIMULATION_STARTED
    time_log = Dict{Symbol, Any}()
    problem = _get_diff_eq_problem(sim.internal, get_simulation_inputs(sim), tspan)
    solution,
    time_log[:timed_solve_time],
    time_log[:solve_bytes_alloc],
    time_log[:sec_in_gc] = @timed SciMLBase.solve(
        problem,
        solver;
        callback = sim.internal.callbacks,
        tstops = sim.internal.tstops,
        progress = get(kwargs, :enable_progress_bar, _prog_meter_enabled()),
        progress_steps = 1,
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
end

"""
    execute!(
        sim::Simulation,
        solver,
        tspan;
        kwargs...
    )

Solves the time-domain dynamic simulation model.

# Arguments
- `sim::Simulation` : Initialized simulation object
- `solver` : Solver used for numerical integration. Must be passed correctly depending on the Type of Simulation Model
- `enable_progress_bar::Bool` : Default: `true`. Enables progress bar for the integration routine.
- Additional solver keyword arguments can be included. See [Common Solver Options](https://diffeq.sciml.ai/stable/basics/common_solver_opts/) in the `DifferentialEquations.jl` documentation for more details.
"""
function execute!(sim::Simulation, solver, tspan; kwargs...)
    logger = configure_logging(sim, "a"; kwargs...)
    try
        Logging.with_logger(logger) do
            _execute!(sim, solver, tspan; kwargs...)
        end
    catch e
        @error "Execution failed" exception = (e, catch_backtrace())
        return sim.status = SIMULATION_FAILED
    finally
        close(logger)
        return sim.status
    end
end

function read_results(sim::Simulation)
    return sim.results
end

function get_dynamic_wrapper(sim::Simulation, name::AbstractString)
    return get_dynamic_wrapper(get_simulation_inputs(sim), name)
end

"""
    get_setpoints(sim::Simulation)

Function that returns the reference setpoints for all the dynamic devices.

# Arguments

- `sim::Simulation` : Simulation object that contains the initial condition and setpoints.
"""
function get_setpoints(sim::Simulation)
    return get_setpoints(get_simulation_inputs(sim))
end
