struct SimulationInputs
    sys::System
    counts::Base.ImmutableDict{Symbol, Int}
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64},Int64}
    dyn_lines::Bool
    voltage_buses_ix::Vector{Int}
    current_buses_ix::Vector{Int}
    global_index::Dict{String,Dict{Symbol,Int64}}
    total_shunts::Dict{Int64,Float64}()
    global_vars::Dict{Symbol,Number}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    aux_arrays::Dict{Int, Vector}
end

function SimulationInputs(;
    sys::PSY.System,
    counts::Base.ImmutableDict{Symbol, Int},
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64},Int64},
    dyn_lines::Bool,
    voltage_buses_ix::Vector{Int} = Vector{Int}(),
    current_buses_ix::Vector{Int} = Vector{Int}(),
    global_index::Dict{String,Dict{Symbol,Int64}} = Dict{String,Dict{Symbol,Int64}}(),
    total_shunts::Dict{Int64,Float64} = Dict{Int64,Float64}(),
    global_vars::Dict{Symbol,Number} = Dict{Symbol,Number}(),
    lookup::Dict{Int, Int} = Dict{Int, Int}(),
    DAE_vector::Vector{Bool} = Vector{Bool}(),
    aux_arrays::Dict{Int, Vector} = Dict{Int, Vector}(),
    )
return SimulationInputs(
    sys,
    counts,
    voltage_buses_ix,
    current_buses_ix,
    global_index,
    Ybus,
    total_shunts,
    global_vars,
    dyn_lines,
    lookup,
    DAE_vector,
    aux_arrays
)
end

get_system(inputs::SimulationInputs) = inputs.sys
get_counts(inputs::SimulationInputs) = inputs.counts
get_voltage_buses_ix(inputs::SimulationInputs) = inputs.voltage_buses_ix
get_current_buses_ix(inputs::SimulationInputs) = inputs.current_buses_ix
get_global_index(inputs::SimulationInputs) = inputs.global_index
get_Ybus(inputs::SimulationInputs) = inputs.Ybus
get_total_shunts(inputs::SimulationInputs) = inputs.total_shunts
get_global_vars(inputs::SimulationInputs) = inputs.global_vars
get_dyn_lines(inputs::SimulationInputs) = inputs.dyn_lines
get_lookup(inputs::SimulationInputs) = inputs.lookup
get_DAE_vector(inputs::SimulationInputs) = inputs.DAE_vector
get_aux_array(inputs::SimulationInputs) = inputs.aux_arrays

get_injection_pointer(inputs::SimulationInputs) =
    get_counts(inputs)[:first_dyn_injection_pointer]
get_branches_pointer(inputs::SimulationInputs) =
    get_counts(inputs)[:first_dyn_branch_point]
get_n_injection_states(inputs::SimulationInputs) = get_counts(inputs)[:injection_n_states]
get_n_branches_states(inputs::SimulationInputs) = get_counts(inputs)[:branches_n_states]
get_system_state_count(inputs::SimulationInputs) = get_counts(inputs)[:total_states]
get_variable_count(inputs::SimulationInputs) = get_counts(inputs)[:total_variables]
get_device_index(inputs::SimulationInputs, device::D) where {D <: PSY.DynamicInjection} =
    get_global_index(inputs)[device.name]
get_bus_count(inputs::SimulationInputs) = get_counts(inputs)[:bus_count]
get_ω_sys(inputs::SimulationInputs) = get_global_vars(inputs)[:ω_sys]
