struct SimulationInputs
    base_power::Float64
    base_frequency::Float64
    dynamic_injectors_data::Vector{DynamicWrapper{<:PSY.DynamicInjection}}
    static_injectors_data::Vector
    static_load_data::Vector
    dynamic_branches::Vector{BranchWrapper}
    injection_n_states::Int
    branches_n_states::Int
    variable_count::Int
    bus_count::Int
    ode_range::UnitRange{Int}
    ybus_rectangular::SparseArrays.SparseMatrixCSC{Float64, Int}
    dyn_lines::Bool
    total_shunts::SparseArrays.SparseMatrixCSC{Float64, Int}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::LinearAlgebra.Diagonal{Float64}
    global_vars_update_pointers::Dict{Int, Int}

    function SimulationInputs(sys::PSY.System, frequency_reference)
        n_buses = get_n_buses(sys)
        Ybus, lookup = _get_ybus(sys)
        wrapped_branches = _wrap_dynamic_branches(sys, lookup)
        has_dyn_lines = !isempty(wrapped_branches)
        branch_state_counts = 2 * length(wrapped_branches)
        injection_start = 2 * n_buses + branch_state_counts + 1
        wrapped_injectors = _wrap_dynamic_injector_data(sys, lookup, injection_start)
        var_count = wrapped_injectors[end].ix_range[end]
        mass_matrix = _make_mass_matrix(wrapped_injectors, var_count, n_buses)
        DAE_vector = _make_DAE_vector(mass_matrix, var_count, n_buses)
        total_shunts = _make_total_shunts(wrapped_branches, n_buses)
        wrapped_loads = _wrap_loads(sys, lookup)
        wrapped_static_injectors = _wrap_static_injectors(sys, lookup)
        sys_f = PSY.get_frequency(sys)
        _adjust_states!(DAE_vector, mass_matrix, total_shunts, n_buses, sys_f)

        global_vars = _make_global_variable_index(
            wrapped_injectors,
            wrapped_static_injectors,
            frequency_reference,
        )

        _make_global_variable_index(
            wrapped_injectors,
            wrapped_static_injectors,
            frequency_reference,
        )

        new(
            PSY.get_base_power(sys),
            sys_f,
            wrapped_injectors,
            wrapped_static_injectors,
            wrapped_loads,
            wrapped_branches,
            var_count - 2 * n_buses - branch_state_counts,
            branch_state_counts,
            var_count,
            n_buses,
            injection_start:var_count,
            Ybus,
            has_dyn_lines,
            total_shunts,
            lookup,
            DAE_vector,
            mass_matrix,
            global_vars,
        )
    end
end

get_base_power(inputs::SimulationInputs) = inputs.base_power
get_base_frequency(inputs::SimulationInputs) = inputs.base_frequency
get_dynamic_injectors_data(inputs::SimulationInputs) = inputs.dynamic_injectors_data
get_dynamic_branches(inputs::SimulationInputs) = inputs.dynamic_branches
get_static_injectiors_data(inputs::SimulationInputs) = inputs.static_injectors_data
get_voltage_buses_ix(inputs::SimulationInputs) = inputs.voltage_buses_ix
get_current_buses_ix(inputs::SimulationInputs) = inputs.current_buses_ix
get_ybus(inputs::SimulationInputs) = inputs.ybus_rectangular
get_total_shunts(inputs::SimulationInputs) = inputs.total_shunts
get_lookup(inputs::SimulationInputs) = inputs.lookup
get_DAE_vector(inputs::SimulationInputs) = inputs.DAE_vector
get_mass_matrix(inputs::SimulationInputs) = inputs.mass_matrix
has_dyn_lines(inputs::SimulationInputs) = inputs.dyn_lines
get_global_vars_update_pointers(inputs::SimulationInputs) =
    inputs.global_vars_update_pointers

get_injection_n_states(inputs::SimulationInputs) = inputs.injection_n_states
get_branches_n_states(inputs::SimulationInputs) = inputs.branches_n_states
get_variable_count(inputs::SimulationInputs) = inputs.variable_count
# TODO: put this in the struct
get_inner_vars_count(inputs::SimulationInputs) =
    inputs.dynamic_injectors_data[end].inner_vars_index[end]
get_ode_range(inputs::SimulationInputs) = inputs.ode_range
get_bus_count(inputs::SimulationInputs) = inputs.bus_count
get_bus_range(inputs::SimulationInputs) = 1:(2 * inputs.bus_count)

"""
SimulationInputs build function for MassMatrixModels
"""
function SimulationInputs(
    ::Type{MassMatrixModel},
    sys::PSY.System,
    frequency_reference = ReferenceBus,
)
    return SimulationInputs(sys, frequency_reference)
end

"""
SimulationInputs build function for ResidualModels
"""
function SimulationInputs(
    ::Type{ResidualModel},
    sys::PSY.System,
    frequency_reference = ReferenceBus,
)
    return SimulationInputs(sys, frequency_reference)
end

function _wrap_dynamic_injector_data(sys::PSY.System, lookup, injection_start::Int)
    injector_data = get_injectors_with_dynamics(sys)
    isempty(injector_data) && error("System doesn't contain any DynamicInjection devices")
    # TODO: Needs a better container that isn't parametrized on an abstract type
    wrapped_injector = Vector(undef, length(injector_data))
    injection_count = 1
    inner_vars_count = 1

    for (ix, device) in enumerate(injector_data)
        @debug PSY.get_name(device)
        dynamic_device = PSY.get_dynamic_injector(device)
        n_states = PSY.get_n_states(dynamic_device)
        ix_range = range(injection_start, length = n_states)
        ode_range = range(injection_count, length = n_states)
        bus_n = PSY.get_number(PSY.get_bus(device))
        bus_ix = lookup[bus_n]
        inner_vars_range = inner_vars_count:get_inner_vars_count(dynamic_device)
        wrapped_injector[ix] =
            DynamicWrapper(device, bus_ix, ix_range, ode_range, inner_vars_range)
        injection_count += n_states
        injection_start += n_states
        inner_vars_count = inner_vars_range[end]
    end
    return wrapped_injector
end

function _wrap_dynamic_branches(sys::PSY.System, lookup::Dict{Int, Int})
    branches_start = 2 * get_n_buses(sys) + 1
    dynamic_branches = get_dynamic_branches(sys)
    wrapped_branches = Vector{BranchWrapper}(undef, length(dynamic_branches))
    if !isempty(wrapped_branches)
        branches_count = 1
        for (ix, br) in enumerate(dynamic_branches)
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = lookup[from_bus_number]
            bus_ix_to = lookup[to_bus_number]
            ix_range = range(branches_start, length = n_states)
            ode_range = range(branches_count, length = n_states)
            branches_count = branches_count + n_states
            wrapped_branches[ix] =
                BranchWrapper(br, bus_ix_from, bus_ix_to, ix_range, ode_range)
            branches_count += n_states
            branches_start += n_states
        end
    else
        @debug("System doesn't contain Dynamic Branches")
    end
    return wrapped_branches
end

function _wrap_static_injectors(sys::PSY.System, lookup::Dict{Int, Int})
    static_injection_data = get_injection_without_dynamics(sys)
    container = Vector{StaticWrapper}(undef, length(static_injection_data))
    isempty(static_injection_data) && return container
    for (ix, ld) in enumerate(static_injection_data)
        bus_n = PSY.get_number(PSY.get_bus(ld))
        bus_ix = lookup[bus_n]
        container[ix] = StaticWrapper(ld, bus_ix)
    end
    return container
end

function _wrap_loads(sys::PSY.System, lookup::Dict{Int, Int})
    # This needs to change if we implement dynamic load models
    static_load_data = PSY.get_components(PSY.ElectricLoad, sys)
    container = Vector(undef, length(static_load_data))
    isempty(static_load_data) && return container
    for (ix, ld) in enumerate(static_load_data)
        bus_n = PSY.get_number(PSY.get_bus(ld))
        bus_ix = lookup[bus_n]
        container[ix] = StaticWrapper(ld, bus_ix)
    end
    return container
end

function _get_ybus(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    dyn_lines = PSY.get_components(PSY.DynamicBranch, sys)
    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus_ = PSY.Ybus(sys)
        ybus = Ybus_[:, :]
        lookup = Ybus_.lookup[1]
        for br in dyn_lines
            ybus_update!(ybus, br, lookup, -1.0)
        end
        # TODO: Improve performance of building the rectangular YBus
        ybus_rectangular = hcat(vcat(real(ybus), -imag(ybus)), vcat(imag(ybus), real(ybus)))
    else
        ybus_rectangular = SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}(
            zeros(2 * n_buses, 2 * n_buses),
        )
        lookup = Dict{Int.Int}()
    end
    return ybus_rectangular, lookup
end

function _static_injection_inputs!(inputs::SimulationInputs, ::Vector{Int}, sys::PSY.System)
    for s in PSY.get_components(PSY.StaticInjection, sys)
        # Temporary condition needs to be removed for other static injections
        PSY.get_dynamic_injector(s) !== nothing && continue
        index_static_injection(inputs, s)
        _attach_control_refs!(s)
    end
    return
end

function _make_mass_matrix(wrapped_injectors, var_count::Int, n_buses::Int)
    mass_matrix = LinearAlgebra.Diagonal(ones(var_count))
    mass_matrix[1:(2 * n_buses), 1:(2 * n_buses)] .= 0.0
    for d in wrapped_injectors
        device_mass_matrix_entries!(mass_matrix, d)
    end
    return mass_matrix
end

function _init_DAE_vector(var_count::Int, n_buses::Int)
    DAE_vector = trues(var_count)
    DAE_vector[1:(n_buses * 2)] .= false
    return DAE_vector
end

function _make_DAE_vector(mass_matrix::AbstractArray, var_count::Int, n_buses::Int)
    DAE_vector = _init_DAE_vector(var_count, n_buses)
    bus_range = 1:(2 * n_buses)
    ode_range = (2 * n_buses + 1):var_count
    mass_buses = @view mass_matrix[bus_range, bus_range]
    DAE_ode = @view DAE_vector[ode_range]
    mass_ode = @view mass_matrix[ode_range, ode_range]
    for i in eachindex(DAE_vector[bus_range])
        IS.@assert_op DAE_vector[bus_range][i] == (mass_buses[i, i] > 0.0)
    end
    for i in eachindex(DAE_ode)
        DAE_ode[i] = (mass_ode[i, i] > 0.0)
    end
    return DAE_vector
end

function _make_total_shunts(wrapped_branches, n_buses::Int)
    shunts = SparseArrays.SparseMatrixCSC{Float64, Int}(zeros(2 * n_buses, 2 * n_buses))
    if isempty(wrapped_branches)
        return shunts
    else
        for br in wrapped_branches
            bus_ix_from = get_bus_ix_from(br)
            bus_ix_to = get_bus_ix_to(br)
            b_from = PSY.get_b(br).from
            b_to = PSY.get_b(br).to
            shunts[bus_ix_from, bus_ix_from + n_buses] += -b_from
            shunts[bus_ix_to, bus_ix_to + n_buses] += -b_to
            shunts[bus_ix_from + n_buses, bus_ix_from] += b_from
            shunts[bus_ix_to + n_buses, bus_ix_to] += b_to
        end
    end
    return shunts
end

function _adjust_states!(
    DAE_vector::BitVector,
    mass_matrix::LinearAlgebra.Diagonal{Float64},
    total_shunts::SparseArrays.SparseMatrixCSC{Float64, Int},
    n_buses::Int,
    sys_f::Float64,
)
    all(iszero.(total_shunts)) && return
    line_constant = 1 / (2.0 * π * sys_f)
    shunts = total_shunts[1:n_buses, (n_buses + 1):end]
    for (ix, val) in enumerate(shunts)
        if val > 0
            mass_matrix[ix, ix] =
                mass_matrix[ix + n_buses, ix + n_buses] = val * line_constant
            DAE_vector[ix] = DAE_vector[ix + n_buses] = true
        end
    end
    return
end

function _make_global_variable_index(
    wrapped_injectors::Vector,
    static_injection_data::Vector,
    frequency_reference::Type{T},
) where {T <: Union{FixedFrequency, ReferenceBus}}
    global_vars_dict = GLOBAL_VARS_IX()
    global_vars_dict[GLOBAL_VAR_SYS_FREQ_INDEX] = get_frequency_reference(
        frequency_reference,
        wrapped_injectors,
        static_injection_data,
    )
    return global_vars_dict
end
