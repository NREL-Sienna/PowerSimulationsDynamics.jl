################ Abstract Structs #######################################
abstract type DynDevice <: PSY.Component end
abstract type DynamicComponent <: PSY.Component end
abstract type GeneratorComponent <: DynamicComponent end
abstract type InverterComponent <: DynamicComponent end
abstract type DynInjection <: DynDevice end
abstract type DynBranch <: DynDevice end

get_bus_number(d::T) where {T <: DynInjection} = d.bus.number
get_bus_number(d::T) where {T <: PSY.Injection} = d.bus.number
total_device_states(d::T) where {T <: DynDevice} = d.n_states

get_ω_ref(d::T) where {T <: DynInjection} = d.ω_ref
get_P_ref(d::T) where {T <: DynInjection} = d.P_ref
get_V_ref(d::T) where {T <: DynInjection} = d.V_ref

################ Base Structs #######################################
mutable struct DynamicSystem
    buses::Vector{PSY.Bus}
    branches::Vector{<:PSY.Branch}
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64}
    Sbase::Float64
    sys_f::Float64
    dyn_injections::Vector{<:DynInjection}
    dyn_branch::Union{Nothing, Vector{<: DynBranch}}
    injections::Vector{<:PSY.Injection}
    global_state_index::Dict{Symbol, Dict{Symbol, Int64}}
    DAE_vector::Vector{Bool}
    counts::Dict{Symbol, Int64}
    indexed::Bool
    ext::Dict{Symbol, Any}
end

function DynamicSystem(Sbase::Float64, sys_f::Float64)
    return DynamicSystem(Vector{PSY.Bus}(),
                        Vector{PSY.Branch}(),
                         SparseMatrixCSC{Complex{Float64}, Int64}(zeros(1,1)),
                         Sbase,
                         sys_f,
                         Vector{DynInjection}(),
                         nothing,
                         Vector{PSY.Injection}(),
                         Dict{Symbol, Dict{Symbol, Int64}}(),
                         Vector{Bool}(),
                         Dict{Symbol, Int64}(),
                         false,
                         Dict{Symbol, Any}())

end

function add_device!(sys::DynamicSystem, device::DynInjection)
    push!(sys.dyn_injections, device)
    sys.indexed = false
    return
end


function add_device!(sys::DynamicSystem, device::DynBranch)
    if isnothing(sys.dyn_branch)
        sys.dyn_branch = Vector{DynBranch}()
    end
    push!(sys.dyn_branch, device)
    sys.indexed = false
    return
end


function add_device!(sys::DynamicSystem, device::PSY.Injection)
    push!(sys.injections, device)
    sys.indexed = false
    return
end

function add_devices!(sys::DynamicSystem, devices)
    for d in devices
        add_device!(sys, d)
    end
    return
end

function add_network!(sys::DynamicSystem, array_bus::Vector{PSY.Bus}, array_branches::Vector{<:PSY.Branch})
    sys.buses = array_bus
    sys.branches = array_branches
    sys.indexed = false
    return
end

function _index_dynamic_system!(sys::DynamicSystem)

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

function DynamicSystem(buses::Vector{PSY.Bus},
                        branches::Vector{Br},
                        dyn_injections::Vector{Di},
                        injections::Vector{I} ,
                        Sbase::Float64,
                        sys_f::Float64,
                        dyn_branch::Union{Nothing, Vector{Dbr}}=nothing) where  {Dbr <: DynBranch,
                                                                                Br <: PSY.Branch,
                                                                                Di <: DynInjection,
                                                                                I <: PSY.Injection}

    dyn_system = DynamicSystem(Sbase, sys_f)

    for dy_i in dyn_injections
        add_device!(dyn_system, dy_i)
    end

    for injection in injections
        add_device!(dyn_system, injection)
    end

    if !isnothing(dyn_branch)
        for dyn_br in dyn_branch
            add_device!(dyn_system, dyn_br)
        end
    end

    add_network!(dyn_system, buses, branches)

    _index_dynamic_system!(dyn_system)

    return dyn_system

end


get_bus_size(sys::DynamicSystem) = length(sys.buses)
get_total_rows(sys::DynamicSystem) = length(sys.DAE_vector)
get_bus_range(sys::DynamicSystem) = 1:2*get_bus_size(sys)
get_injection_pointer(sys::DynamicSystem) = sys.counts[:first_dyn_injection_pointer]
get_branches_pointer(sys::DynamicSystem) = sys.counts[:first_dyn_branch_point]
get_n_injection_states(sys::DynamicSystem) = sys.counts[:injection_n_states]
get_n_branches_states(sys::DynamicSystem) = sys.counts[:branches_n_states]
system_state_count(sys::DynamicSystem) = sys.counts[:total_states]
get_sys_base(sys::DynamicSystem) = sys.Sbase
get_device_index(sys::DynamicSystem, device::D) where {D <: DynDevice} = sys.global_state_index[device.name]
get_sys_f(sys::DynamicSystem) = sys.sys_f
is_indexed(sys::DynamicSystem) = sys.indexed
