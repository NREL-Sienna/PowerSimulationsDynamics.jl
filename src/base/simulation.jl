    mutable struct Simulation
    system::PSY.System
    problem::DiffEqBase.DAEProblem
    perturbations::Vector{<:Perturbation}
    x0_init::Vector{Float64}
    initialized::Bool
    tstops::Vector{Float64}
    callbacks::Union{Nothing, DiffEqBase.CallbackSet}
    solution::Union{Nothing, DiffEqBase.DAESolution}
    ext::Dict{String, Any}
end

function Simulation(system::PSY.System,
                    tspan::NTuple{2, Float64},
                    perturbations::Vector{<:Perturbation};
                    initialize_simulation::Bool=true,
                    kwargs...)

    initialized = false
    n_buses = length(PSY.get_components(PSY.Bus, system))
    DAE_vector = collect(falses(n_buses*2))
    _index_dynamic_system!(DAE_vector, system)
    var_count = get_variable_count(system)

    if initialize_simulation
        @info("Initializing Simulation States")
        x0_guess = get(kwargs, :initial_guess, zeros(var_count))
        x0_init, initialized = _calculate_initial_conditions(system, x0_guess)
    else
        x0_init = zeros(var_count)
    end

    dx0 = zeros(var_count)
    prob = DiffEqBase.DAEProblem(system_model!,
                                 dx0,
                                 x0_init,
                                 tspan,
                                 (control, system),
                                 differential_vars = DAE_vector)

    callback_set, tstops = _build_callbacks(perturbations::Vector{<:Perturbation})

    return Simulation(system,
                     prob,
                     perturbations,
                     x0_init,
                     initialized,
                     tstops,
                     callback_set,
                     nothing
                     Dict{String, Any}()
                     )

end

function _build_callbacks(perturbations::Vector{<:Perturbation})
    return nothing, [0.0]
end


function _calculate_initial_conditions(sys::PSY.System, initial_guess::Vector{Float64})
    # TODO: Code to refine initial_guess
    var_count = get_variable_count(sys)
    dx0 = zeros(var_count) #Define a vector of zeros for the derivative
    inif! = (out,x) -> system_model!(
        out, #output of the function
        dx0, #derivatives equal to zero
        x, #states
        ([0.0], sys), #Parameters: [0.0] is not used
        0.0) #time equals to zero.
    sys_solve = NLsolve.nlsolve(inif!, initial_guess) #Solve using initial guess x0
    if !NLsolve.converged(sys_solve)
        @warn("Initialization failed, initial conditions do not meet conditions for an stable equilibrium")
    end

    return sys_solve.zero, NLsolve.converged(sys_solve)
end


function run_simulation!(sim::Simulation, solver; kwargs...)
    sim.solution = DiffEqBase.solve(sim.problem, solver;
                           callback = sim.callbacks,
                           tstops = sim.tstops,
                           kwargs...)
    return
end

function _index_local_states!(component_state_index::Vector{Int64},
                              local_states::Vector{Symbol},
                              component::PSY.DynamicComponent)
    for (ix, s) in enumerate(component.states)
        component_state_index[ix] = findfirst(x->x == s, local_states)
    end
    return
end

function _attach_ports!(component::PSY.DynamicComponent)
    component.ext[PORTS] = Ports(component)
    return
end

function _attach_inner_vars!(device::PSY.DynamicGenerator)
    device.ext[INNER_VARS] = zeros(8)
    return
end

function _attach_inner_vars!(device::PSY.DynamicInverter)
    device.ext[INNER_VARS] = zeros(13)
    return
end

function _index_port_mapping!(index_component_inputs::Vector{Int64},
                             local_states::Vector{Symbol},
                              component::PSY.DynamicComponent)
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

function _index_dynamic_system!(DAE_vector::Vector{Bool},
                                sys::PSY.System)

    global_state_index = Dict{String, Dict{Symbol, Int64}}()
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    state_space_ix = n_buses*2
    total_states = 0
    first_dyn_branch_point = -1
    branches_n_states = 0

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        if !(:states in fieldnames(typeof(d)))
            continue
        end
        _make_device_index!(d)
        device_n_states = PSY.get_n_states(d)
        DAE_vector = vcat(DAE_vector, collect(trues(device_n_states)))
        total_states += device_n_states
        state_ix = Dict{Symbol, Int}()
        for s in PSY.get_states(d)
            state_space_ix += 1
            state_ix[s] = state_space_ix
        end
        global_state_index[PSY.get_name(d)] = state_ix
    end
    injection_n_states = state_space_ix - n_buses*2

#=
    if !(isnothing(sys.dyn_branch))
        first_dyn_branch_point =  state_space_ix + 1
        for br in sys.dyn_branch
            arc = br.arc
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            sys.DAE_vector[from_bus_number] = sys.DAE_vector[from_bus_number+n_buses] = true
            sys.DAE_vector[to_bus_number] = sys.DAE_vector[to_bus_number+n_buses] = true
            sys.DAE_vector = vcat(sys.DAE_vector, collect(trues(br.n_states)))
            total_states += br.n_states
            state_ix = Dict{Symbol, Int}()
            for (ix, s) in enumerate(getfield(br, :states))
                state_space_ix += 1
                state_ix[s] = state_space_ix
            end
            sys.global_state_index[getfield(br, :name)] = state_ix
        end


        for (ix, val) in enumerate(sys.DAE_vector[1:n_buses])
            if val
                sys.global_state_index[Symbol("V_$(ix)")] = Dict(:R => ix,
                                                                :I => ix + n_buses)
                total_states += 2
                state_space_ix += 2
            end
        end
        branches_n_states = state_space_ix - injection_n_states - n_buses*2
    end
=#
    @assert total_states == state_space_ix - n_buses*2

    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus = PSY.Ybus(sys)[:, :]
    else
        Ybus = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(n_buses, n_buses))
    end
    sys_ext = Dict{String, Any}() #I change it to be Any
    counts = Dict{Symbol, Int64}(:total_states => total_states,
                                  :injection_n_states => injection_n_states,
                                  :branches_n_states => branches_n_states,
                                  :first_dyn_injection_pointer => 2*n_buses+1,
                                  :first_dyn_branch_point => first_dyn_branch_point,
                                  :total_variables => total_states + 2*n_buses)

    sys_ext[LITS_COUNTS] = counts
    sys_ext[GLOBAL_INDEX] = global_state_index
    sys_ext[YBUS] = Ybus
    sys.internal.ext = sys_ext

    return
end

get_injection_pointer(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:first_dyn_injection_pointer]
get_branches_pointer(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:first_dyn_branch_point]
get_n_injection_states(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:injection_n_states]
get_n_branches_states(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:branches_n_states]
get_system_state_count(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:total_states]
get_variable_count(sys::PSY.System) = PSY.get_ext(sys)[LITS_COUNTS][:total_variables]
get_device_index(
    sys::PSY.System,
    device::D,
    ) where {D <: PSY.DynamicInjection} = PSY.get_ext(sys)[GLOBAL_INDEX][device.name]

get_inner_vars(device::PSY.DynamicInjection) = device.ext[INNER_VARS]

function _get_internal_mapping(device::PSY.DynamicInjection, key::AbstractString, ty::Type{T}) where T <: PSY.DynamicComponent
    device_index = PSY.get_ext(device)[key]
    val = get(device_index, ty, nothing)
    @assert !isnothing(val)
    return val
end

function get_local_state_ix(device::PSY.DynamicInjection, ty::Type{T}) where T <: PSY.DynamicComponent
    return _get_internal_mapping(device, LOCAL_STATE_MAPPING, ty)
end


function get_input_port_ix(device::PSY.DynamicInjection, ty::Type{T}) where T <: PSY.DynamicComponent
    return _get_internal_mapping(device, INPUT_PORT_MAPPING, ty)
end
