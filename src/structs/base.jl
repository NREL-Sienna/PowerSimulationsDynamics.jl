
################ Abstract Structs #######################################
abstract type DynDevice <: PSY.Component end
abstract type DynamicComponent <: PSY.Component end
abstract type GeneratorComponent <: DynamicComponent end
abstract type InverterComponent <: DynamicComponent end
abstract type DynInjection <: DynDevice end
abstract type DynBranch <: DynDevice end

get_bus_number(d::T) where {T <: DynInjection} = d.bus.number
get_bus_number(d::T) where {T <: PSY.Injection} = d.bus.number #
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
    ext::Dict{Symbol, Any}
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
                n_buses = length(buses)
                total_states = 0
                state_space_ix = 0
                global_state_index = Dict{Symbol, Dict{Symbol, Int64}}()
                DAE_vector = collect(falses(n_buses*2))
                total_states = 0
                state_space_ix = n_buses*2
                first_dyn_branch_point = -1
                branches_n_states = 0

                for d in dyn_injections
                    if !(:states in fieldnames(typeof(d)))
                        continue
                    end
                    device_n_states = getfield(d, :n_states)
                    DAE_vector = vcat(DAE_vector, collect(trues(device_n_states)))
                    total_states += device_n_states
                    state_ix = Dict{Symbol, Int}()
                    for s in getfield(d, :states)
                        state_space_ix += 1
                        state_ix[s] = state_space_ix
                    end
                    global_state_index[getfield(d, :name)] = state_ix
                end
                injection_n_states = state_space_ix - n_buses*2

                if !(isnothing(dyn_branch))
                    first_dyn_branch_point =  state_space_ix + 1
                    for br in dyn_branch
                        arc = br.arc
                        from_bus_number = PSY.get_number(arc.from)
                        to_bus_number = PSY.get_number(arc.to)
                        DAE_vector[from_bus_number] = DAE_vector[from_bus_number+n_buses] = true
                        DAE_vector[to_bus_number] = DAE_vector[to_bus_number+n_buses] = true
                        DAE_vector = vcat(DAE_vector, collect(trues(br.n_states)))
                        total_states += br.n_states
                        state_ix = Dict{Symbol, Int}()
                        for (ix, s) in enumerate(getfield(br, :states))
                            state_space_ix += 1
                            state_ix[s] = state_space_ix
                        end
                        global_state_index[getfield(br, :name)] = state_ix
                    end


                    for (ix, val) in enumerate(DAE_vector[1:n_buses])
                        if val
                            global_state_index[Symbol("V_$(ix)")] = Dict(:R => ix,
                                                                         :I => ix + n_buses)
                            total_states += 2
                            state_space_ix += 2
                        end
                    end
                    branches_n_states = state_space_ix - injection_n_states - n_buses*2
                end

            @assert total_states == state_space_ix - n_buses*2

            if !isempty(branches)
                Ybus = PSY.Ybus(branches, buses)[:, :]
            else
                Ybus = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(n_buses, n_buses))
            end

            counts = Dict{Symbol, Int64}(:total_states => total_states,
                                         :injection_n_states => injection_n_states,
                                         :branches_n_states => branches_n_states,
                                         :first_dyn_injection_pointer => 2*n_buses+1,
                                         :first_dyn_branch_point => first_dyn_branch_point)
            ext = Dict{Symbol, Any}()
            new(buses,
                branches,
                Ybus,
                Sbase,
                sys_f,
                dyn_injections,
                dyn_branch,
                injections,
                global_state_index,
                DAE_vector,
                counts,
                ext)
        end
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
