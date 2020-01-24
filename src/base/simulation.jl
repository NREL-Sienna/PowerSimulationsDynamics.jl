mutable struct Simulation
    system::PSY.System
    problem::DiffEqBase.DAEProblem
    #callbacks::DiffEqBase.DiscreteCallback
    #tstops::Vector{Float64}
    x0_init::Vector{Float64}
    initialized::Bool
    solution::Union{Nothing, DiffEqBase.DAESolution}
    #global_state_index::Dict{Symbol, Dict{Symbol, Int64}}
    #DAE_vector::Vector{Bool}
    #counts::Dict{Symbol, Int64}
    #indexed::Bool
    #ext::Dict{Symbol, Any}
end

function Simulation(dyn_system::PSY.System,
                           tspan,
                           #control,
                           #callback,
                           #x0_init
                           )

    if !is_indexed(dyn_system)
        _index_dynamic_system!(dyn_system)
    end

    dx0 = zeros(get_total_rows(dyn_system))
    prob = DiffEqBase.DAEProblem(system_model,
                              dx0,
                              x0_init,
                              tspan,
                              (control, dyn_system),
                              differential_vars = dyn_system.DAE_vector)

    return Simulation(dyn_system,
                             prob,
                             callback,
                             [1.0],
                             x0_init,
                             true,
                             nothing)

end


function run_simulation!(sim::Simulation, solver; kwargs...)
    sim.solution = DiffEqBase.solve(sim.problem,
                           solver;
                           callback = sim.callbacks,
                           tstops = sim.tstops, kwargs...)

    return
end

function _make_local_state_mapping!(device::PSY.DynamicInjection)

    local_state_mapping = Dict{PSY.DynamicComponent, Vector{Int64}}()
    local_state_space_ix = 0
    for c in PSY.get_dynamic_components(device)
        component_state_index= Vector{Int64}(undef, length(c.states))
        for (ix, s) in enumerate(c.states)
            component_state_index[ix] = findfirst(x->x == s, states)
        end
        local_state_mapping[c] = component_state_index
    end

    device.ext["local_state_mapping"] = local_state_mapping

    return
end

function _make_port_mapping(device::PSY.DynamicInjection)
    states = PSY.get_states(device)
    input_port_mapping = Dict{PSY.DynamicComponent, Vector{Int64}}()

    for c in components
        index_component_in = Vector{Int64}()
        for i in c.ports.state
            tmp = [(ix,var) for (ix,var) in enumerate(states) if var == i]
            isempty(tmp) && continue
            push!(index_component_in, tmp[1][1])
        end
        input_port_mapping[c] = index_component_in
    end

    return input_port_mapping

end

function _index_dynamic_system!(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    total_states = 0
    state_space_ix = 0
    DAE_vector = collect(falses(n_buses*2))
    total_states = 0
    state_space_ix = n_buses*2
    first_dyn_branch_point = -1
    branches_n_states = 0
    global_state_index = Dict{String, Dict{Symbol, Int64}}()

    for d in PSY.get_components(PSY.DynamicInjection, sys)
        if !(:states in fieldnames(typeof(d)))
            continue
        end
        _make_local_state_mapping!(d)
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

    counts = Dict{Symbol, Int64}(:total_states => total_states,
                                    :injection_n_states => injection_n_states,
                                    :branches_n_states => branches_n_states,
                                    :first_dyn_injection_pointer => 2*n_buses+1,
                                    :first_dyn_branch_point => first_dyn_branch_point)

    indexed = true

    return

end
