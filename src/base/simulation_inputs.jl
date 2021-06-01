#TODO: Make inmmutable later. Requires refactor of ThreePhase Fault callbacks
mutable struct SimulationInputs
    sys::PSY.System
    injectors_data::Vector{<:PSY.StaticInjection}
    dynamic_branches::Vector{PSY.DynamicBranch}
    counts::Base.ImmutableDict{Symbol, Int}
    Ybus::SparseMatrixCSC{Complex{Float64}, Int}
    dyn_lines::Bool
    voltage_buses_ix::Vector{Int}
    current_buses_ix::Vector{Int}
    global_index::Dict{String, Dict{Symbol, Int}}
    total_shunts::Dict{Int, Float64}
    global_vars::Dict{Symbol, Number}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::SparseMatrixCSC{Float64, Int}
    aux_arrays::Dict{Int, Vector}
    tspan::NTuple{2, Float64}
end

function SimulationInputs(;
    sys::PSY.System,
    counts::Base.ImmutableDict{Symbol, Int} = Base.ImmutableDict{Symbol, Int}(),
    voltage_buses_ix::Vector{Int} = Vector{Int}(),
    global_index::Dict{String, Dict{Symbol, Int}} = MAPPING_DICT(),
    total_shunts::Dict{Int, Float64} = Dict{Int, Float64}(),
    global_vars::Dict{Symbol, Number} = GLOBAL_VARS_IX(),
    DAE_vector::Vector{Bool} = Vector{Bool}(),
    mass_matrix::SparseMatrixCSC{Float64, Int} = SparseMatrixCSC{Float64, Int}(zeros(1, 1)),
    aux_arrays::Dict{Int, Vector} = Dict{Int, Vector}(),
    tspan::NTuple{2, Float64} = (0.0, 0.0),
)
    Ybus, lookup = _get_Ybus(sys)
    n_buses = length(PSY.get_bus_numbers(sys))

    injector_data = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing && PSY.get_available(x),
    )

    if isempty(injector_data)
        error("System doesn't contain any DynamicInjection devices")
    end

    dynamic_branches =
        collect(PSY.get_components(PSY.DynamicBranch, sys, x -> PSY.get_available(x)))

    return SimulationInputs(
        sys,
        collect(injector_data),
        dynamic_branches,
        counts,
        Ybus,
        !isempty(dynamic_branches),
        voltage_buses_ix,
        collect(1:n_buses),
        global_index,
        total_shunts,
        global_vars,
        lookup,
        DAE_vector,
        mass_matrix,
        aux_arrays,
        tspan,
    )
end


function add_aux_arrays!(inputs::SimulationInputs, ::Type{T}) where {T <: Number}
    @debug "Auxiliary Arrays created with Type $(T)"
    bus_count = get_bus_count(inputs)
    get_aux_arrays(inputs)[1] = collect(zeros(T, bus_count))                       #I_injections_r
    get_aux_arrays(inputs)[2] = collect(zeros(T, bus_count))                       #I_injections_i
    get_aux_arrays(inputs)[3] = collect(zeros(T, get_n_injection_states(inputs)))  #injection_ode
    get_aux_arrays(inputs)[4] = collect(zeros(T, get_n_branches_states(inputs)))   #branches_ode
    get_aux_arrays(inputs)[5] = collect(zeros(Complex{T}, bus_count))              #I_bus
    get_aux_arrays(inputs)[6] = collect(zeros(T, 2 * bus_count))                   #I_balance
    return
end

function _init_DAE_vector!(inputs::SimulationInputs, n_buses::Int)
    push!(inputs.DAE_vector, collect(falses(n_buses * 2))...)
end

function _dynamic_lines_inputs!(
    inputs::SimulationInputs,
    state_space_ix::Vector{Int},
    n_buses::Int,
)
    static_bus_var_count = 2 * n_buses
    branches_n_states = 0
    if inputs.dyn_lines
        global_state_index = inputs.global_index
        for br in get_dynamic_branches(inputs)
            index_dynamic_lines!(inputs, br, n_buses)
            add_states_to_global!(global_state_index, state_space_ix, br)
            branches_n_states += PSY.get_n_states(br)
        end

        for (ix, val) in enumerate(inputs.DAE_vector[1:n_buses])
            if val
                #This V_ix should be V_number.
                global_state_index["V_$(ix)"] = Dict(:R => ix, :I => ix + n_buses)
                static_bus_var_count -= 2
                push!(inputs.voltage_buses_ix, ix)
                @assert static_bus_var_count >= 0
            end
        end
        unique!(inputs.voltage_buses_ix)
        setdiff!(inputs.current_buses_ix, inputs.voltage_buses_ix)
    else
        @debug("System doesn't contain Dynamic Branches")
    end
    return branches_n_states, static_bus_var_count
end

function _static_injection_inputs!(
    inputs::SimulationInputs,
    ::Vector{Int},
    sys::PSY.System,
)
    for s in PSY.get_components(PSY.StaticInjection, sys)
        index_static_injection(inputs, s)
    end
    return
end

function _dynamic_injection_inputs!(
    inputs::SimulationInputs,
    state_space_ix::Vector{Int},
)
    dynamic_injection_states = 0
    for d in get_injectors_data(inputs)
        @debug PSY.get_name(d)
        _attach_control_refs!(d)
        dynamic_injector = PSY.get_dynamic_injector(d)
        dynamic_injection_states += PSY.get_n_states(dynamic_injector)
        index_dynamic_injection(inputs, dynamic_injector, state_space_ix)
    end
    return dynamic_injection_states
end

function build!(inputs::SimulationInputs, ::Type{ImplicitModel})
    sys = get_system(inputs)
    n_buses = length(PSY.get_bus_numbers(sys))
    _init_DAE_vector!(inputs, n_buses)
    state_space_ix = Int[n_buses * 2]
    branches_n_states, static_bus_var_count =
        _dynamic_lines_inputs!(inputs, state_space_ix, n_buses)

    _static_injection_inputs!(inputs, state_space_ix, sys)
    injection_n_states = _dynamic_injection_inputs!(inputs, state_space_ix)

    set_frequency_reference!(inputs, sys)
    IS.@assert_op get_ω_sys(inputs) != -1

    IS.@assert_op injection_n_states == state_space_ix[1] - branches_n_states - n_buses * 2
    IS.@assert_op n_buses * 2 - static_bus_var_count >= 0

    inputs.counts = Base.ImmutableDict(
        :total_states => state_space_ix[1] - static_bus_var_count,
        :injection_n_states => injection_n_states,
        :branches_n_states => branches_n_states,
        :first_dyn_injection_pointer => 2 * n_buses + branches_n_states + 1,
        :first_dyn_branch_point => 2 * n_buses + 1,
        :total_variables => state_space_ix[1],
        :bus_count => n_buses,
    )

    @debug inputs.counts

    return inputs
end

get_system(inputs::SimulationInputs) = inputs.sys
get_injectors_data(inputs::SimulationInputs) = inputs.injectors_data
get_dynamic_branches(inputs::SimulationInputs) = inputs.dynamic_branches
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
get_aux_arrays(inputs::SimulationInputs) = inputs.aux_arrays
get_tspan(inputs::SimulationInputs) = inputs.tspan

get_injection_pointer(inputs::SimulationInputs) =
    get_counts(inputs)[:first_dyn_injection_pointer]
get_branches_pointer(inputs::SimulationInputs) = get_counts(inputs)[:first_dyn_branch_point]
get_n_injection_states(inputs::SimulationInputs) = get_counts(inputs)[:injection_n_states]
get_n_branches_states(inputs::SimulationInputs) = get_counts(inputs)[:branches_n_states]
get_system_state_count(inputs::SimulationInputs) = get_counts(inputs)[:total_states]
get_variable_count(inputs::SimulationInputs) = get_counts(inputs)[:total_variables]
get_device_index(inputs::SimulationInputs, device::D) where {D <: PSY.DynamicInjection} =
    get_global_index(inputs)[device.name]
get_bus_count(inputs::SimulationInputs) = get_counts(inputs)[:bus_count]
get_ω_sys(inputs::SimulationInputs) = get_global_vars(inputs)[:ω_sys]
