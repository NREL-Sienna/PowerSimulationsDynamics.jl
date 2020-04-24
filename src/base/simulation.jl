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
    ext::Dict{String, Any}
end

function Simulation(
    system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbations::Vector{<:Perturbation} = Vector{Perturbation}();
    initialize_simulation::Bool = true,
    kwargs...,
)
    check_kwargs(kwargs, SIMULATION_ACCEPTED_KWARGS, "Simulation")
    initialized = false
    DAE_vector = _index_dynamic_system!(system)
    var_count = get_variable_count(system)

    flat_start = zeros(var_count)
    bus_count = length(PSY.get_components(PSY.Bus, system))
    flat_start[1:bus_count] .= 1.0
    x0_init = get(kwargs, :initial_guess, flat_start)

    if initialize_simulation
        @info("Initializing Simulation States")
        initialized = _calculate_initial_conditions!(system, x0_init)
    end

    dx0 = zeros(var_count)
    callback_set, tstops = _build_perturbations(perturbations::Vector{<:Perturbation})

    _add_aux_arrays!(system, Float64)
    prob = DiffEqBase.DAEProblem(
        system!,
        dx0,
        x0_init,
        tspan,
        system,
        differential_vars = DAE_vector;
        kwargs...,
    )

    return Simulation(
        system,
        false,
        prob,
        perturbations,
        x0_init,
        initialized,
        tstops,
        callback_set,
        nothing,
        Dict{String, Any}(),
    )
end

function Simulation(
    system::PSY.System,
    tspan::NTuple{2, Float64},
    perturbation::Perturbation;
    initialize_simulation::Bool = true,
    kwargs...,
)
    return Simulation(
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
        6 => collect(zeros(T, 2 * bus_count)),                      #I_balance
    )
    system.internal.ext[AUX_ARRAYS] = aux_arrays
    return
end

function _build_perturbations(perturbations::Vector{<:Perturbation})
    isempty(perturbations) && return DiffEqBase.CallbackSet(), [0.0]
    perturbations_count = length(perturbations)
    callback_vector = Vector{DiffEqBase.DiscreteCallback}(undef, perturbations_count)
    tstops = Vector{Float64}(undef, perturbations_count)
    for (ix, pert) in enumerate(perturbations)
        condition = (x, t, integrator) -> t in [pert.time]
        affect = get_affect(pert)
        callback_vector[ix] = DiffEqBase.DiscreteCallback(condition, affect)
        tstops[ix] = pert.time
    end
    callback_tuple = Tuple(cb for cb in callback_vector)
    callback_set = DiffEqBase.CallbackSet((), callback_tuple)
    return callback_set, tstops
end

function _calculate_initial_conditions!(sys::PSY.System, initial_guess::Vector{Float64})
    # TODO: Code to refine initial_guess
    _add_aux_arrays!(sys, Real)
    bus_count = length(PSY.get_components(PSY.Bus, sys))
    var_count = get_variable_count(sys)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    inif! = (out, x) -> system!(
        out,    #output of the function
        dx0,    #derivatives equal to zero
        x,      #states
        sys,    #Parameters
        0.0,    #time equals to zero.
    )
    sys_solve = NLsolve.nlsolve(
        inif!,
        initial_guess,
        xtol = :1e-9,
        ftol = :1e-9,
        method = :trust_region,
    ) #Solve using initial guess x0
    if !NLsolve.converged(sys_solve)
        @warn("Initialization failed, initial conditions do not meet conditions for an stable equilibrium")
    end
    initial_guess .= sys_solve.zero
    return NLsolve.converged(sys_solve)
end

function _index_local_states!(
    component_state_index::Vector{Int64},
    local_states::Vector{Symbol},
    component::PSY.DynamicComponent,
)
    for (ix, s) in enumerate(component.states)
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
    device.ext[INNER_VARS] = zeros(T, 13)
    return
end

function _attach_control_refs!(device::PSY.DynamicInjection)
    device.ext[CONTROL_REFS] = [
        PSY.get_V_ref(device),
        PSY.get_ω_ref(device),
        PSY.get_P_ref(device),
        PSY.get_Q_ref(device),
    ]
    return
end

function _index_port_mapping!(
    index_component_inputs::Vector{Int64},
    local_states::Vector{Symbol},
    component::PSY.DynamicComponent,
)
    _attach_ports!(component)
    for i in component.ext[PORTS].state
        tmp = [(ix, var) for (ix, var) in enumerate(local_states) if var == i]
        isempty(tmp) && continue
        push!(index_component_inputs, tmp[1][1])
    end

    return
end

function _make_device_index!(device::PSY.DynamicInjection)
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
    global_state_index::Dict{String, Dict{Symbol, Int64}},
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
    global_state_index = Dict{String, Dict{Symbol, Int64}}()
    state_space_ix = [n_buses * 2]
    current_buses_no = collect(1:n_buses)
    static_bus_var_count = 2 * length(current_buses_no)
    voltage_buses_no = Vector{Int}()
    total_states = 0
    first_dyn_branch_point = -1
    branches_n_states = 0
    global_vars = Dict{Symbol, Number}(
        :ω_sys => 1.0,
        :ω_sys_index => -1, #To define 0 if infinite source, bus_number otherwise,
    )
    total_shunts = Dict{Int, Float64}()
    found_ref_bus = false

    dyn_branches = PSY.get_components(DynamicLine, sys)
    if !(isempty(dyn_branches))
        first_dyn_branch_point = state_space_ix[1] + 1
        for br in dyn_branches
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            merge!(
                +,
                total_shunts,
                Dict(
                    from_bus_number => 1 / PSY.get_b(br).from,
                    to_bus_number => 1 / PSY.get_b(br).to,
                ),
            )
            push!(voltage_buses_no, from_bus_number, to_bus_number)
            DAE_vector[from_bus_number] = DAE_vector[from_bus_number + n_buses] = true
            DAE_vector[to_bus_number] = DAE_vector[to_bus_number + n_buses] = true
            DAE_vector = push!(DAE_vector, collect(trues(n_states))...)
            total_states += n_states
            _add_states_to_global!(global_state_index, state_space_ix, br)
        end

        for (ix, val) in enumerate(DAE_vector[1:n_buses])
            if val
                global_state_index["V_$(ix)"] = Dict(:R => ix, :I => ix + n_buses)
                total_states += 2
                static_bus_var_count -= 2
                push!(voltage_buses_no, ix)
                @assert static_bus_var_count >= 0
            end
        end
        branches_n_states = state_space_ix[1] - n_buses * 2
    else
        @debug("System doesn't contain Dynamic Branches")
    end
    unique!(voltage_buses_no)
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
        _make_device_index!(d)
        device_n_states = PSY.get_n_states(d)
        DAE_vector = push!(DAE_vector, collect(trues(device_n_states))...)
        total_states += device_n_states
        _add_states_to_global!(global_state_index, state_space_ix, d)
        btype != PSY.BusTypes.REF && continue
        global_vars[:ω_sys_index] = global_state_index[d.name][:ω] #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end
    injection_n_states = state_space_ix[1] - branches_n_states - n_buses * 2
    @assert total_states == state_space_ix[1] - static_bus_var_count
    @debug total_states
    setdiff!(current_buses_no, voltage_buses_no)
    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus = PSY.Ybus(sys)[:, :]
    else
        Ybus = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(n_buses, n_buses))
    end
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

    sys_ext[LITS_COUNTS] = counts
    sys_ext[GLOBAL_INDEX] = global_state_index
    sys_ext[VOLTAGE_BUSES_NO] = voltage_buses_no
    sys_ext[CURRENT_BUSES_NO] = current_buses_no
    sys_ext[YBUS] = Ybus
    sys_ext[TOTAL_SHUNTS] = total_shunts
    sys_ext[GLOBAL_VARS] = global_vars
    @assert sys_ext[GLOBAL_VARS][:ω_sys_index] != -1
    sys.internal.ext = sys_ext
    return DAE_vector
end

get_injection_pointer(sys::PSY.System) =
    PSY.get_ext(sys)[LITS_COUNTS][:first_dyn_injection_pointer]
get_branches_pointer(sys::PSY.System) =
    PSY.get_ext(sys)[LITS_COUNTS][:first_dyn_branch_point]
get_n_injection_states(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:injection_n_states]
get_n_branches_states(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:branches_n_states]
get_system_state_count(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:total_states]
get_variable_count(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:total_variables]
get_device_index(sys::PSY.System, device::D) where {D <: PSY.DynamicInjection} =
    PSY.get_ext(sys)[GLOBAL_INDEX][device.name]
get_inner_vars(device::PSY.DynamicInjection) = device.ext[INNER_VARS]
get_ω_sys(sys::PSY.System) = PSY.get_ext(sys)[GLOBAL_VARS][:ω_sys]
get_current_bus_no(sys::PSY.System) = PSY.get_ext(sys)[CURRENT_BUSES_NO]
get_voltage_bus_no(sys::PSY.System) = PSY.get_ext(sys)[VOLTAGE_BUSES_NO]
get_total_shunts(sys::PSY.System) = PSY.get_ext(sys)[TOTAL_SHUNTS]
get_bus_count(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:bus_count]

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
    for real_eig in real(vals)
        real_eig >= 0.0 && return false
    end
    return true
end

function small_signal_analysis(sim::Simulation; kwargs...)
    if sim.reset
        @error("Reset the simulation")
    end
    _change_vector_type(sim.system)
    _add_aux_arrays!(sim.system, Real)
    var_count = LITS.get_variable_count(sim.system)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    bus_count = length(PSY.get_components(PSY.Bus, sim.system))
    sysf! = (out, x) -> LITS.system!(
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
    for b_no in get_voltage_bus_no(sim.system)
        diff_states[b_no] = true
        diff_states[b_no + n_buses] = true
    end
    alg_states = .!diff_states
    fx = @view jacobian[diff_states, diff_states]
    gy = jacobian[alg_states, alg_states]
    fy = @view jacobian[diff_states, alg_states]
    gx = @view jacobian[alg_states, diff_states]
    # TODO: Make operation using BLAS!
    reduced_jacobian = fx - fy * inv(gy) * gx
    vals, vect = LinearAlgebra.eigen(reduced_jacobian)
    return SmallSignalOutput(
        reduced_jacobian,
        vals,
        vect,
        _determine_stability(vals),
        x_eval,
    )
end
