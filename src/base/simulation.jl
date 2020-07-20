mutable struct Simulation
    system::PSY.System
    reset::Bool
    problem::DiffEqBase.DAEProblem
    perturbations::Vector{<:Perturbation}
    x0_init::Vector{Float64}
    initialized::Bool
    tstops::Vector{Float64}
    callbacks::DiffEqBase.CallbackSet
    solution::Union{Nothing, DiffEqBase.DAESolution}
    simulation_folder::String
    ext::Dict{String, Any}
end

"""
Initializes the simulations and builds the indexing. The initial conditions are stored in the system.

# Accepted Key Words
- `system_to_file::Bool`: Serializes the initialized system
"""
function Simulation!(
    simulation_folder::String,
    system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
)
    check_folder(simulation_folder)
    sim = build_simulation(simulation_folder, system, tspan, perturbations; kwargs...)
    if get(kwargs, :system_to_file, false)
        PSY.to_json(system, joinpath(simulation_folder, "initialized_system.json"))
    end
    return sim
end

"""
Initializes the simulations and builds the indexing. The input system is not modified during the initialization

# Accepted Key Words
- `system_to_file::Bool`: Serializes the original input system
"""
function Simulation(
    simulation_folder::String,
    system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
)
    check_folder(simulation_folder)
    simulation_system = deepcopy(system)
    if get(kwargs, :system_to_file, false)
        PSY.to_json(system, joinpath(simulation_folder, "input_system.json"),)
    end
    return build_simulation(simulation_folder, simulation_system, tspan, perturbations; kwargs...)
end

function build_simulation(simulation_folder::String,
    simulation_system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    kwargs...,
)
    PSY.set_units_base_system!(simulation_system, "device_base")
    check_kwargs(kwargs, SIMULATION_ACCEPTED_KWARGS, "Simulation")
    initialized = false
    DAE_vector = _index_dynamic_system!(simulation_system)
    var_count = get_variable_count(simulation_system)

    flat_start = zeros(var_count)
    bus_count = length(PSY.get_components(PSY.Bus, simulation_system))
    flat_start[1:bus_count] .= 1.0
    x0_init = get(kwargs, :initial_guess, flat_start)

    initialize_simulation = get(kwargs, :initialize_simulation, true)
    if initialize_simulation
        @info("Initializing Simulation States")
        _add_aux_arrays!(simulation_system, Real)
        initialized = calculate_initial_conditions!(simulation_system, x0_init)
    end

    dx0 = zeros(var_count)
    callback_set, tstops = _build_perturbations(simulation_system, perturbations)
    _add_aux_arrays!(simulation_system, Float64)
    prob = DiffEqBase.DAEProblem(
        system!,
        dx0,
        x0_init,
        tspan,
        simulation_system,
        differential_vars = DAE_vector;
        kwargs...,
    )
    return Simulation(
        simulation_system,
        false,
        prob,
        perturbations,
        x0_init,
        initialized,
        tstops,
        callback_set,
        nothing,
        simulation_folder,
        Dict{String, Any}(),
    )
end

function Simulation(
    simulation_folder::String,
    system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbation::Perturbation;
    initialize_simulation::Bool = true,
    kwargs...,
)
    return Simulation(
        simulation_folder,
        system,
        tspan,
        [perturbation];
        initialize_simulation = initialize_simulation,
        kwargs...,
    )
end

function _add_aux_arrays!(system::PSY.System, T)
    bus_count = get_bus_count(system)
    aux_arrays = Dict(
        1 => collect(zeros(T, bus_count)),                       #I_injections_r
        2 => collect(zeros(T, bus_count)),                       #I_injections_i
        3 => collect(zeros(T, get_n_injection_states(system))),  #injection_ode
        4 => collect(zeros(T, get_n_branches_states(system))),   #branches_ode
        5 => collect(zeros(Complex{T}, bus_count)),              #I_bus
        6 => collect(zeros(T, 2 * bus_count)),                   #I_balance
    )
    system.internal.ext[AUX_ARRAYS] = aux_arrays
    return
end

function _build_perturbations(system::PSY.System, perturbations::Vector{<:Perturbation})
    isempty(perturbations) && return DiffEqBase.CallbackSet(), [0.0]
    perturbations_count = length(perturbations)
    callback_vector = Vector{DiffEqBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Vector{Float64}(undef, perturbations_count)
    for (ix, pert) in enumerate(perturbations)
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(system, pert)
        callback_vector[ix] = DiffEqBase.DiscreteCallback(condition, affect)
        tstops[ix] = pert.time
    end
    callback_set = DiffEqBase.CallbackSet((), tuple(callback_vector...))
    return callback_set, tstops
end

function _index_local_states!(
    component_state_index::Vector{Int64},
    local_states::Vector{Symbol},
    component::PSY.DynamicComponent,
)
    for (ix, s) in enumerate(PSY.get_states(component))
        component_state_index[ix] = findfirst(x -> x == s, local_states)
    end
    return
end

function _attach_ports!(component::PSY.DynamicComponent)
    component.ext[PORTS] = Ports(component)
    return
end

function _attach_inner_vars!(
    device::PSY.DynamicGenerator,
    ::Type{T} = Real,
) where {T <: Real}
    device.ext[INNER_VARS] = zeros(T, 8)
    return
end

function _attach_inner_vars!(
    device::PSY.DynamicInverter,
    ::Type{T} = Real,
) where {T <: Real}
    device.ext[INNER_VARS] = zeros(T, 14)
    return
end

function _attach_control_refs!(device::PSY.DynamicInjection)
    device.ext[CONTROL_REFS] = [
        PSY.get_V_ref(device),
        PSY.get_ω_ref(device),
        PSY.get_P_ref(device),
        PSY.get_reactive_power(PSY.get_static_injector(device)),
    ]
    return
end

function _index_port_mapping!(
    index_component_inputs::Vector{Int64},
    local_states::Vector{Symbol},
    component::PSY.DynamicComponent,
)
    _attach_ports!(component)
    for i in component.ext[PORTS].states
        tmp = [(ix, var) for (ix, var) in enumerate(local_states) if var == i]
        isempty(tmp) && continue
        push!(index_component_inputs, tmp[1][1])
    end

    return
end

function _get_Ybus(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    dyn_lines = PSY.get_components(PSY.DynamicBranch, sys)
    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus_ = PSY.Ybus(sys)
        Ybus = Ybus_[:, :]
        lookup = Ybus_.lookup[1]
        for br in dyn_lines
            ybus!(Ybus, br, lookup, -1.0)
        end
    else
        Ybus = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(n_buses, n_buses))
        lookup = Dict{Int.Int}()
    end
    return Ybus, lookup
end

function _make_device_index!(device::PSY.DynamicInjection, sys_basepower::Float64)
    states = PSY.get_states(device)
    device_state_mapping = Dict{Type{<:PSY.DynamicComponent}, Vector{Int64}}()
    input_port_mapping = Dict{Type{<:PSY.DynamicComponent}, Vector{Int64}}()
    _attach_inner_vars!(device)
    _attach_control_refs!(device)

    for c in PSY.get_dynamic_components(device)
        device_state_mapping[typeof(c)] = Vector{Int64}(undef, length(c.states))
        input_port_mapping[typeof(c)] = Vector{Int64}()
        _index_local_states!(device_state_mapping[typeof(c)], states, c)
        _index_port_mapping!(input_port_mapping[typeof(c)], states, c)
        device.ext[LOCAL_STATE_MAPPING] = device_state_mapping
        device.ext[INPUT_PORT_MAPPING] = input_port_mapping
    end

    return
end

function _add_states_to_global!(
    global_state_index::MAPPING_DICT,
    state_space_ix::Vector{Int64},
    device::PSY.Device,
)
    device_state_ix = Dict{Symbol, Int}()
    for s in PSY.get_states(device)
        state_space_ix[1] += 1
        device_state_ix[s] = state_space_ix[1]
    end
    global_state_index[PSY.get_name(device)] = device_state_ix
    return
end

function _index_dynamic_system!(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    DAE_vector = collect(falses(n_buses * 2))
    global_state_index = MAPPING_DICT()
    state_space_ix = [n_buses * 2]
    current_buses_ix = collect(1:n_buses)
    static_bus_var_count = 2 * length(current_buses_ix)
    voltage_buses_ix = Vector{Int}()
    total_states = 0
    first_dyn_branch_point = -1
    branches_n_states = 0
    global_vars = Dict{Symbol, Number}(
        :ω_sys => 1.0,
        :ω_sys_index => -1, #To define 0 if infinite source, bus_number otherwise,
    )
    total_shunts = Dict{Int, Float64}()
    found_ref_bus = false
    sys_basepower = PSY.get_base_power(sys)

    Ybus, lookup = _get_Ybus(sys)
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sys)

    if !(isempty(dyn_branches))
        first_dyn_branch_point = state_space_ix[1] + 1
        for br in dyn_branches
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = lookup[from_bus_number]
            bus_ix_to = lookup[to_bus_number]
            merge!(
                +,
                total_shunts,
                Dict(
                    bus_ix_from => 1 / PSY.get_b(br).from,
                    bus_ix_to => 1 / PSY.get_b(br).to,
                ),
            )
            push!(voltage_buses_ix, bus_ix_from, bus_ix_to)
            DAE_vector[bus_ix_from] = DAE_vector[bus_ix_from + n_buses] = true
            DAE_vector[bus_ix_to] = DAE_vector[bus_ix_to + n_buses] = true
            DAE_vector = push!(DAE_vector, collect(trues(n_states))...)
            total_states += n_states
            _add_states_to_global!(global_state_index, state_space_ix, br)
        end

        for (ix, val) in enumerate(DAE_vector[1:n_buses])
            if val
                #This V_ix should be V_number.
                global_state_index["V_$(ix)"] = Dict(:R => ix, :I => ix + n_buses)
                total_states += 2
                static_bus_var_count -= 2
                push!(voltage_buses_ix, ix)
                @assert static_bus_var_count >= 0
            end
        end
        branches_n_states = state_space_ix[1] - n_buses * 2
    else
        @debug("System doesn't contain Dynamic Branches")
    end
    unique!(voltage_buses_ix)
    sources = PSY.get_components(PSY.Source, sys)
    for s in sources
        btype = PSY.get_bustype(PSY.get_bus(s))
        if (btype == PSY.BusTypes.REF) && found_ref_bus
            throw(IS.ConflictingInputsError("The system can't have more than one source or generator in the REF Bus"))
        end
        btype != PSY.BusTypes.REF && continue
        global_vars[:ω_sys_index] = 0 #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end
    dynamic_injection = PSY.get_components(PSY.DynamicInjection, sys)
    isempty(dynamic_injection) &&
        error("System doesn't contain any DynamicInjection devices")
    for d in dynamic_injection
        if !(:states in fieldnames(typeof(d)))
            continue
        end
        device_bus = PSY.get_bus(d)
        btype = PSY.get_bustype(device_bus)
        if (btype == PSY.BusTypes.REF) && found_ref_bus
            throw(IS.ConflictingInputsError("The system can't have more than one source or generator in the REF Bus"))
        end
        _make_device_index!(d, sys_basepower)
        device_n_states = PSY.get_n_states(d)
        DAE_vector = push!(DAE_vector, collect(trues(device_n_states))...)
        total_states += device_n_states
        _add_states_to_global!(global_state_index, state_space_ix, d)
        btype != PSY.BusTypes.REF && continue
        global_vars[:ω_sys_index] = global_state_index[PSY.get_name(d)][:ω] #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end
    injection_n_states = state_space_ix[1] - branches_n_states - n_buses * 2
    @assert total_states == state_space_ix[1] - static_bus_var_count
    @debug total_states
    setdiff!(current_buses_ix, voltage_buses_ix)
    sys_ext = Dict{String, Any}()
    counts = Base.ImmutableDict(
        :total_states => total_states,
        :injection_n_states => injection_n_states,
        :branches_n_states => branches_n_states,
        :first_dyn_injection_pointer => 2 * n_buses + branches_n_states + 1,
        :first_dyn_branch_point => first_dyn_branch_point,
        :total_variables => total_states + static_bus_var_count,
        :bus_count => n_buses,
    )

    sys_ext[PSID_COUNTS] = counts
    sys_ext[GLOBAL_INDEX] = global_state_index
    sys_ext[VOLTAGE_BUSES_IX] = voltage_buses_ix
    sys_ext[CURRENT_BUSES_IX] = current_buses_ix
    sys_ext[YBUS] = Ybus
    sys_ext[TOTAL_SHUNTS] = total_shunts
    sys_ext[GLOBAL_VARS] = global_vars
    sys_ext[DYN_LINES] = !isempty(PSY.get_components(PSY.DynamicBranch, sys))
    sys_ext[LOOKUP] = lookup
    @assert sys_ext[GLOBAL_VARS][:ω_sys_index] != -1
    sys.internal.ext = sys_ext
    return DAE_vector
end

get_injection_pointer(sys::PSY.System) =
    PSY.get_ext(sys)[PSID_COUNTS][:first_dyn_injection_pointer]
get_branches_pointer(sys::PSY.System) =
    PSY.get_ext(sys)[PSID_COUNTS][:first_dyn_branch_point]
get_n_injection_states(sys::PSY.System) = PSY.get_ext(sys)[PSID_COUNTS][:injection_n_states]
get_n_branches_states(sys::PSY.System) = PSY.get_ext(sys)[PSID_COUNTS][:branches_n_states]
get_system_state_count(sys::PSY.System) = PSY.get_ext(sys)[PSID_COUNTS][:total_states]
get_variable_count(sys::PSY.System) = PSY.get_ext(sys)[PSID_COUNTS][:total_variables]
get_device_index(sys::PSY.System, device::D) where {D <: PSY.DynamicInjection} =
    PSY.get_ext(sys)[GLOBAL_INDEX][device.name]
get_inner_vars(device::PSY.DynamicInjection) = device.ext[INNER_VARS]
get_ω_sys(sys::PSY.System) = PSY.get_ext(sys)[GLOBAL_VARS][:ω_sys]
get_current_bus_ix(sys::PSY.System) = PSY.get_ext(sys)[CURRENT_BUSES_IX]
get_voltage_bus_ix(sys::PSY.System) = PSY.get_ext(sys)[VOLTAGE_BUSES_IX]
get_total_shunts(sys::PSY.System) = PSY.get_ext(sys)[TOTAL_SHUNTS]
get_bus_count(sys::PSY.System) = PSY.get_ext(sys)[PSID_COUNTS][:bus_count]
get_lookup(sys::PSY.System) = PSY.get_ext(sys)[LOOKUP]

function _get_internal_mapping(
    device::PSY.DynamicInjection,
    key::AbstractString,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    device_index = PSY.get_ext(device)[key]
    val = get(device_index, ty, nothing)
    @assert !isnothing(val)
    return val
end

function get_local_state_ix(
    device::PSY.DynamicInjection,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    return _get_internal_mapping(device, LOCAL_STATE_MAPPING, ty)
end

function get_input_port_ix(
    device::PSY.DynamicInjection,
    ty::Type{T},
) where {T <: PSY.DynamicComponent}
    return _get_internal_mapping(device, INPUT_PORT_MAPPING, ty)
end

function run_simulation!(sim::Simulation, solver; kwargs...)
    if sim.reset
        @error("Reset the simulation")
    end

    sim.solution = DiffEqBase.solve(
        sim.problem,
        solver;
        callback = sim.callbacks,
        tstops = sim.tstops,
        kwargs...,
    )
    return
end

function _change_vector_type(sys::PSY.System)
    for d in PSY.get_components(PSY.DynamicInjection, sys)
        _attach_inner_vars!(d, Real)
    end
end

function _determine_stability(vals::Vector{Complex{Float64}})
    stable = true
    for real_eig in real(vals)
        real_eig > 0.0 && return false
    end
    return true
end

function small_signal_analysis(sim::Simulation; kwargs...)
    if sim.reset
        @error("Reset the simulation")
    end
    _change_vector_type(sim.system)
    _add_aux_arrays!(sim.system, Real)
    var_count = PSID.get_variable_count(sim.system)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    bus_count = length(PSY.get_components(PSY.Bus, sim.system))
    sysf! = (out, x) -> system!(
        out,            #output of the function
        dx0,            #derivatives equal to zero
        x,              #states
        sim.system,     #Parameters
        0.0,            #time equals to zero.
    )
    out = zeros(var_count) #Define a vector of zeros for the output
    x_eval = get(kwargs, :operating_point, sim.x0_init)
    jacobian = ForwardDiff.jacobian(sysf!, out, x_eval)
    n_buses = length(PSY.get_components(PSY.Bus, sim.system))
    diff_states = collect(trues(var_count))
    diff_states[1:(2 * n_buses)] .= false
    for b_ix in get_voltage_bus_ix(sim.system)
        diff_states[b_ix] = true
        diff_states[b_ix + n_buses] = true
    end
    alg_states = .!diff_states
    fx = @view jacobian[diff_states, diff_states]
    gy = jacobian[alg_states, alg_states]
    fy = @view jacobian[diff_states, alg_states]
    gx = @view jacobian[alg_states, diff_states]
    # TODO: Make operation using BLAS!
    reduced_jacobian = fx - fy * inv(gy) * gx
    vals, vect = LinearAlgebra.eigen(reduced_jacobian)
    sources = collect(PSY.get_components(PSY.Source, sim.system))
    if isempty(sources)
        @warn("No Infinite Bus found. Confirm stability directly checking eigenvalues.\nIf all eigenvalues are on the left-half plane and only one eigenvalue is zero, the system is small signal stable.")
        info_evals = "Eigenvalues are:\n"
        for i in vals
            info_evals = info_evals * string(i) * "\n"
        end
        @info(info_evals)
    end
    return SmallSignalOutput(
        reduced_jacobian,
        vals,
        vect,
        _determine_stability(vals),
        x_eval,
    )
end
