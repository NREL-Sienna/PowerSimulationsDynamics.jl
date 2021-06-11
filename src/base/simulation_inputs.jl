struct SimulationInputs
    base_power::Float64
    base_frequency::Float64
    injectors_data::Vector{<:PSY.StaticInjection}
    static_injection_data::Vector{<:PSY.StaticInjection}
    dynamic_branches::Vector{PSY.DynamicBranch}
    counts::Base.ImmutableDict{Symbol, Int}
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}
    dyn_lines::Bool
    voltage_buses_ix::Vector{Int}
    current_buses_ix::Vector{Int}
    global_index::Dict{String, Dict{Symbol, Int}}
    total_shunts::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}
    global_vars::Dict{Symbol, Real}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::SparseArrays.SparseMatrixCSC{Float64, Int}
    aux_arrays::Dict{Int, Vector}
    tspan::NTuple{2, Float64}
end

function SimulationInputs(sys::PSY.System, tspan::NTuple{2, Float64} = (0.0, 0.0))
    injector_data = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing && PSY.get_available(x),
    )

    if isempty(injector_data)
        error("System doesn't contain any DynamicInjection devices")
    end

    Ybus, lookup = _get_Ybus(sys)
    n_buses = length(PSY.get_bus_numbers(sys))

    dynamic_branches =
        collect(PSY.get_components(PSY.DynamicBranch, sys, x -> PSY.get_available(x)))

    branch_state_counts = 2 * length(dynamic_branches)
    injector_state_count = sum(PSY.get_n_states.(PSY.get_dynamic_injector.(injector_data)))
    var_count = injector_state_count + 2 * n_buses + branch_state_counts

    mass_matrix = SparseArrays.sparse(LinearAlgebra.I, var_count, var_count)
    mass_matrix[1:(2 * n_buses), 1:(2 * n_buses)] .= 0.0

    static_injection_data = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) === nothing && PSY.get_available(x),
    )

    counts = Base.ImmutableDict(
        :injection_n_states => injector_state_count,
        :branches_n_states => branch_state_counts,
        :first_dyn_injection_pointer => 2 * n_buses + branch_state_counts + 1,
        :first_dyn_branch_point => isempty(dynamic_branches) ? 0 : 2 * n_buses + 1,
        :total_variables => var_count,
        :bus_count => n_buses,
    )

    return SimulationInputs(
        PSY.get_base_power(sys),
        PSY.get_frequency(sys),
        collect(injector_data),
        collect(static_injection_data),
        dynamic_branches,
        counts,
        Ybus,
        !isempty(dynamic_branches),
        Vector{Int}(),
        collect(1:n_buses),
        MAPPING_DICT(),
        SparseArrays.spzeros(Complex{Float64}, n_buses, n_buses),
        GLOBAL_VARS_IX(),
        lookup,
        _init_DAE_vector!(var_count, n_buses),
        mass_matrix,
        Dict{Int, Vector{Real}}(),
        tspan,
    )
end

get_base_power(inputs::SimulationInputs) = inputs.base_power
get_base_frequency(inputs::SimulationInputs) = inputs.base_frequency
get_injectors_data(inputs::SimulationInputs) = inputs.injectors_data
get_dynamic_branches(inputs::SimulationInputs) = inputs.dynamic_branches
get_static_injections_data(inputs::SimulationInputs) = inputs.static_injection_data
get_counts(inputs::SimulationInputs) = inputs.counts
get_voltage_buses_ix(inputs::SimulationInputs) = inputs.voltage_buses_ix
get_current_buses_ix(inputs::SimulationInputs) = inputs.current_buses_ix
get_global_index(inputs::SimulationInputs) = inputs.global_index
get_Ybus(inputs::SimulationInputs) = inputs.Ybus
get_total_shunts(inputs::SimulationInputs) = inputs.total_shunts
get_global_vars(inputs::SimulationInputs) = inputs.global_vars
get_lookup(inputs::SimulationInputs) = inputs.lookup
get_DAE_vector(inputs::SimulationInputs) = inputs.DAE_vector
get_mass_matrix(inputs::SimulationInputs) = inputs.mass_matrix
get_aux_arrays(inputs::SimulationInputs) = inputs.aux_arrays
get_tspan(inputs::SimulationInputs) = inputs.tspan
has_dyn_lines(inputs::SimulationInputs) = inputs.dyn_lines

get_injection_pointer(inputs::SimulationInputs) =
    get_counts(inputs)[:first_dyn_injection_pointer]
get_branches_pointer(inputs::SimulationInputs) = get_counts(inputs)[:first_dyn_branch_point]
get_n_injection_states(inputs::SimulationInputs) = get_counts(inputs)[:injection_n_states]
get_n_branches_states(inputs::SimulationInputs) = get_counts(inputs)[:branches_n_states]
get_variable_count(inputs::SimulationInputs) = get_counts(inputs)[:total_variables]
get_device_index(inputs::SimulationInputs, device::D) where {D <: PSY.DynamicInjection} =
    get_global_index(inputs)[device.name]
get_bus_count(inputs::SimulationInputs) = get_counts(inputs)[:bus_count]
get_ω_sys(inputs::SimulationInputs) = get_global_vars(inputs)[:ω_sys]

function change_vector_type!(inputs::SimulationInputs, ::Type{T}) where {T <: Real}
    for d in get_injectors_data(inputs)
        attach_inner_vars!(PSY.get_dynamic_injector(d), T)
    end
    return
end

function simulation_pre_step!(inputs::SimulationInputs, ::Type{T}) where {T <: Real}
    add_aux_arrays!(inputs, T)
    change_vector_type!(inputs, T)
    return
end

function add_aux_arrays!(inputs::SimulationInputs, ::Type{T}) where {T <: Real}
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

function _get_Ybus(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    dyn_lines = PSY.get_components(PSY.DynamicBranch, sys)
    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus_ = PSY.Ybus(sys)
        Ybus = Ybus_[:, :]
        lookup = Ybus_.lookup[1]
        for br in dyn_lines
            ybus_update!(Ybus, br, lookup, -1.0)
        end
    else
        Ybus = SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}(zeros(n_buses, n_buses))
        lookup = Dict{Int.Int}()
    end
    return Ybus, lookup
end

function _init_DAE_vector!(n_vars::Int, n_buses::Int)
    DAE_vector = trues(n_vars)
    DAE_vector[1:(n_buses * 2)] .= false
    return DAE_vector
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

function _static_injection_inputs!(inputs::SimulationInputs, ::Vector{Int}, sys::PSY.System)
    for s in PSY.get_components(PSY.StaticInjection, sys)
        index_static_injection(inputs, s)
    end
    return
end

function _dynamic_injection_inputs!(inputs::SimulationInputs, state_space_ix::Vector{Int})
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

function _mass_matrix_inputs!(inputs::SimulationInputs)
    mass_matrix = get_mass_matrix(inputs)
    injection_range = PSID.get_injection_pointer(inputs):PSID.get_variable_count(inputs)
    mass_matrix_injectors = view(mass_matrix, injection_range, injection_range)
    for d in PSID.get_injectors_data(inputs)
        dynamic_injector = PSY.get_dynamic_injector(d)
        device_mass_matrix_entries!(mass_matrix_injectors, dynamic_injector)
    end

    if has_dyn_lines(inputs)
        sys_f = get_base_frequency(inputs)
        shunts = get_total_shunts(inputs)
        n_buses = get_bus_count(inputs)
        for i in 1:n_buses
            val = imag(shunts[i, i])
            if val > 0
                mass_matrix[i, i] =
                    mass_matrix[i + n_buses, i + n_buses] = val * (1 / (2.0 * π * sys_f))
            end
        end
    end
end

# Default implementation for both models. This implementation is to future proof if there is
# a divergence between the required build methods
function _build!(inputs::SimulationInputs, sys::PSY.System)
    n_buses = get_bus_count(inputs)
    state_space_ix = Int[n_buses * 2]
    branches_n_states, static_bus_var_count =
        _dynamic_lines_inputs!(inputs, state_space_ix, n_buses)

    _static_injection_inputs!(inputs, state_space_ix, sys)
    injection_n_states = _dynamic_injection_inputs!(inputs, state_space_ix)

    set_frequency_reference!(inputs, sys)
    IS.@assert_op get_ω_sys(inputs) != -1

    _mass_matrix_inputs!(inputs)

    IS.@assert_op injection_n_states == state_space_ix[1] - branches_n_states - n_buses * 2
    IS.@assert_op n_buses * 2 - static_bus_var_count >= 0
    IS.@assert_op length(inputs.DAE_vector) == state_space_ix[1]

    @debug inputs.counts
    return
end

"""
SimulationInputs build function for MassMatrixModels
"""
function build!(inputs::SimulationInputs, ::Type{MassMatrixModel}, sys::PSY.System)
    return _build!(inputs, sys)
end

"""
SimulationInputs build function for ImplicitModels
"""
function build!(inputs::SimulationInputs, ::Type{ImplicitModel}, sys::PSY.System)
    return _build!(inputs, sys)
end
