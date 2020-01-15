mutable struct DynamicSimulation
    dyn_system::PSY.System
    problem::DiffEqBase.DAEProblem
    callbacks::DiffEqBase.DiscreteCallback
    tstops::Vector{Float64}
    x0_init::Vector{Float64}
    initialized::Bool
    solution::Union{Nothing, DiffEqBase.DAESolution}
    #global_state_index::Dict{Symbol, Dict{Symbol, Int64}}
    #DAE_vector::Vector{Bool}
    #counts::Dict{Symbol, Int64}
    #indexed::Bool
    #ext::Dict{Symbol, Any}
end

function DynamicSimulation(dyn_system::PSY.System,
                           tspan,
                           control,
                           callback,
                           x0_init)

    if !is_indexed(dyn_system)
        _index_dynamic_system!(dyn_system)
    end

    dx0 = zeros(get_total_rows(dyn_system))
    prob = DiffEqBase.DAEProblem(system_model!,
                              dx0,
                              x0_init,
                              tspan,
                              (control, dyn_system),
                              differential_vars = dyn_system.DAE_vector)

    return DynamicSimulation(dyn_system,
                             prob,
                             callback,
                             [1.0],
                             x0_init,
                             true,
                             nothing)

end


function run_simulation!(sim::DynamicSimulation, solver; kwargs...)
    sim.solution = DiffEqBase.solve(sim.problem,
                           solver;
                           callback = sim.callbacks,
                           tstops = sim.tstops, kwargs...)

    return
end


function _index_dynamic_system!(sys::DynamicSimulation)

    n_buses = length(sys.buses)
    total_states = 0
    state_space_ix = 0
    sys.DAE_vector = collect(falses(n_buses*2))
    total_states = 0
    state_space_ix = n_buses*2
    first_dyn_branch_point = -1
    branches_n_states = 0

    for d in sys.dyn_injections
        if !(:states in fieldnames(typeof(d)))
            continue
        end
        device_n_states = getfield(d, :n_states)
        sys.DAE_vector = vcat(sys.DAE_vector, collect(trues(device_n_states)))
        total_states += device_n_states
        state_ix = Dict{Symbol, Int}()
        for s in getfield(d, :states)
            state_space_ix += 1
            state_ix[s] = state_space_ix
        end
        sys.global_state_index[getfield(d, :name)] = state_ix
    end
    injection_n_states = state_space_ix - n_buses*2

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

    @assert total_states == state_space_ix - n_buses*2

    if !isempty(sys.branches)
        sys.Ybus = PSY.Ybus(sys.branches, sys.buses)[:, :]
    else
        sys.Ybus = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(n_buses, n_buses))
    end

    sys.counts = Dict{Symbol, Int64}(:total_states => total_states,
                                    :injection_n_states => injection_n_states,
                                    :branches_n_states => branches_n_states,
                                    :first_dyn_injection_pointer => 2*n_buses+1,
                                    :first_dyn_branch_point => first_dyn_branch_point)

    sys.indexed = true

    return

end
