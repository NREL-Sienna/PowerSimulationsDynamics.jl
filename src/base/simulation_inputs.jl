struct SimulationInputs
    dynamic_injectors::Vector{DynamicWrapper{<:PSY.DynamicInjection}}
    static_injectors::Vector
    static_loads::Vector
    dynamic_branches::Vector{BranchWrapper}
    injection_n_states::Int
    branches_n_states::Int
    variable_count::Int
    inner_vars_count::Int
    bus_count::Int
    parameter_count::Int
    ode_range::UnitRange{Int}
    dyn_lines::Bool
    lookup::Dict{Int, Int}
    global_vars_update_pointers::Dict{Int, Int}
    global_state_map::MAPPING_DICT
    global_inner_var_map::Dict{String, Dict}
    ybus_rectangular::SparseArrays.SparseMatrixCSC{Float64, Int}
    total_shunts::SparseArrays.SparseMatrixCSC{Float64, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::LinearAlgebra.Diagonal{Float64}
    parameters::Vector{Float64}
    constant_lags::Vector
end

function SimulationInputs(
    sys::PSY.System,
    ::T,
) where {T <: Union{ConstantFrequency, ReferenceBus}}
    n_buses = get_n_buses(sys)
    _, lookup = _get_ybus(sys)
    state_count = 2 * n_buses + 1
    parameter_count = 1
    TimerOutputs.@timeit BUILD_TIMER "Wrap Branches" begin
        wrapped_branches, state_count, parameter_count =
            _wrap_dynamic_branches(sys, lookup, state_count, parameter_count)
        has_dyn_lines = !isempty(wrapped_branches)
    end
    n_branch_states = state_count - (2 * n_buses + 1)
    injection_start = state_count

    TimerOutputs.@timeit BUILD_TIMER "Wrap Dynamic Injectors" begin
        wrapped_injectors, state_count, parameter_count =
            _wrap_dynamic_injector_data(sys, lookup, state_count, parameter_count)
        n_vars = state_count - 1
    end

    TimerOutputs.@timeit BUILD_TIMER "Wrap Static Injectors" begin
        wrapped_loads, parameter_count = _wrap_loads(sys, lookup, parameter_count)
        wrapped_static_injectors, parameter_count =
            _wrap_static_injectors(sys, lookup, parameter_count)
        n_parameters = parameter_count - 1
    end
    global_vars =
        _make_global_variable_index(wrapped_injectors, wrapped_static_injectors, T)

    inner_vars_count = 0
    for i in length(wrapped_injectors):-1:1
        if length(wrapped_injectors[i].inner_vars_index) > 0
            inner_vars_count = wrapped_injectors[i].inner_vars_index[end]
            break
        end
    end
    Ybus, _ = _get_ybus(sys)
    total_shunts = _make_total_shunts(wrapped_branches, n_buses)
    TimerOutputs.@timeit BUILD_TIMER "Build initial parameters" begin
        parameter_count = n_parameters
        initial_parameters = zeros(parameter_count)
        _update_initial_parameters!(initial_parameters, wrapped_branches)
        _update_initial_parameters!(initial_parameters, wrapped_injectors)
        _update_initial_parameters!(initial_parameters, wrapped_loads)
        _update_initial_parameters!(initial_parameters, wrapped_static_injectors)
    end
    mass_matrix = _make_mass_matrix(wrapped_injectors, n_vars, n_buses)
    DAE_vector = _make_DAE_vector(mass_matrix, n_vars, n_buses)
    _adjust_states!(
        DAE_vector,
        mass_matrix,
        total_shunts,
        n_buses,
        PSY.get_frequency(sys),
    )

    delays = get_constant_lags(sys)
    if !isempty(delays)
        @info "System has delays. Use the correct solver for delay differential equations."
    end

    return SimulationInputs(
        wrapped_injectors,
        wrapped_static_injectors,
        wrapped_loads,
        wrapped_branches,
        n_vars - 2 * n_buses - n_branch_states,
        n_branch_states,
        n_vars,
        inner_vars_count,
        n_buses,
        n_parameters,
        injection_start:n_vars,
        has_dyn_lines,
        lookup,
        global_vars,
        MAPPING_DICT(),
        Dict{String, Dict}(),
        Ybus,
        total_shunts,
        DAE_vector,
        mass_matrix,
        initial_parameters,
        delays,
    )
end
#LEVEL ONE INPUTS
get_dynamic_injectors(inputs::SimulationInputs) = inputs.dynamic_injectors
get_dynamic_branches(inputs::SimulationInputs) = inputs.dynamic_branches
get_static_injectors(inputs::SimulationInputs) = inputs.static_injectors
get_static_loads(inputs::SimulationInputs) = inputs.static_loads
get_lookup(inputs::SimulationInputs) = inputs.lookup
has_dyn_lines(inputs::SimulationInputs) = inputs.dyn_lines
get_global_vars_update_pointers(inputs::SimulationInputs) =
    inputs.global_vars_update_pointers
get_global_state_map(inputs::SimulationInputs) = inputs.global_state_map
get_injection_n_states(inputs::SimulationInputs) = inputs.injection_n_states
get_branches_n_states(inputs::SimulationInputs) = inputs.branches_n_states
get_variable_count(inputs::SimulationInputs) = inputs.variable_count
get_inner_vars_count(inputs::SimulationInputs) = inputs.inner_vars_count
get_ode_ouput_range(inputs::SimulationInputs) = inputs.ode_range
get_bus_count(inputs::SimulationInputs) = inputs.bus_count
get_parameter_count(inputs::SimulationInputs) = inputs.parameter_count
get_bus_range(inputs::SimulationInputs) = 1:(2 * inputs.bus_count)

#LEVEL TWO INPUTS
get_ybus(inputs::SimulationInputs) = inputs.ybus_rectangular
get_total_shunts(inputs::SimulationInputs) = inputs.total_shunts

#LEVEL THREE INPUTS
get_DAE_vector(inputs::SimulationInputs) = inputs.DAE_vector
get_mass_matrix(inputs::SimulationInputs) = inputs.mass_matrix
get_parameters(inputs::SimulationInputs) = inputs.parameters
get_constant_lags(inputs::SimulationInputs) = inputs.constant_lags

# Utility function not to be used for performance sensitive operations
function get_voltage_buses_ix(inputs::SimulationInputs)
    n_buses = get_bus_count(inputs::SimulationInputs)
    return findall(get_DAE_vector(inputs)[1:n_buses])
end

# Utility function not to be used for performance sensitive operations
function get_current_buses_ix(inputs::SimulationInputs)
    n_buses = get_bus_count(inputs::SimulationInputs)
    return findall(.!get_DAE_vector(inputs)[1:n_buses])
end

"""
SimulationInputs build function for MassMatrixModels
"""
function SimulationInputs(
    ::Type{MassMatrixModel},
    sys::PSY.System,
    frequency_reference::Union{ConstantFrequency, ReferenceBus},
)
    return SimulationInputs(sys, frequency_reference)
end

"""
SimulationInputs build function for ResidualModels
"""
function SimulationInputs(
    ::Type{ResidualModel},
    sys::PSY.System,
    frequency_reference::Union{ConstantFrequency, ReferenceBus},
)
    return SimulationInputs(sys, frequency_reference)
end

function _update_initial_parameters!(initial_parameters, wrapped_devices)
    for wrapped_device in wrapped_devices
        p = get_params(wrapped_device)
        p_range = get_p_range(wrapped_device)
        initial_parameters[p_range] .= p
    end
end

function _get_n_params(
    dynamic_device::PSY.DynamicInjection,
    static_device::PSY.StaticInjection,
)
    return get_n_params(dynamic_device)
end

function _get_n_params(
    dynamic_device::T,
    static_device::PSY.StaticInjection,
) where {T <: Union{PSY.DynamicGenerator, PSY.DynamicInverter}}
    return get_n_params(dynamic_device) + get_n_params(static_device)
end

function _wrap_dynamic_injector_data(
    sys::PSY.System,
    lookup,
    state_count::Int,
    parameter_count::Int,
)   #injection_start
    injector_data = get_injectors_with_dynamics(sys)
    isempty(injector_data) && error("System doesn't contain any DynamicInjection devices")
    # TODO: Needs a better container that isn't parametrized on an abstract type
    wrapped_injector = Vector(undef, length(injector_data))
    injection_count = 1
    inner_vars_count = 1
    sys_base_power = PSY.get_base_power(sys)
    sys_base_freq = PSY.get_frequency(sys)
    @assert !isempty(injector_data)
    for (ix, device) in enumerate(injector_data)
        @debug "Wrapping $(PSY.get_name(device))"
        dynamic_device = PSY.get_dynamic_injector(device)
        n_states = PSY.get_n_states(dynamic_device)
        n_params = _get_n_params(dynamic_device, device)
        n_inner_vars = get_inner_vars_count(dynamic_device)
        ix_range = range(state_count; length = n_states)
        ode_range = range(injection_count; length = n_states)
        p_range = range(parameter_count; length = n_params)
        bus_n = PSY.get_number(PSY.get_bus(device))
        bus_ix = lookup[bus_n]
        inner_vars_range = range(inner_vars_count; length = n_inner_vars)
        @debug "ix_range=$ix_range ode_range=$ode_range inner_vars_range= $inner_vars_range p_range=$p_range"
        dynamic_device = PSY.get_dynamic_injector(device)
        @assert dynamic_device !== nothing
        wrapped_injector[ix] = DynamicWrapper(
            device,
            dynamic_device,
            bus_ix,
            ix_range,
            ode_range,
            p_range,
            inner_vars_range,
            sys_base_power,
            sys_base_freq,
        )
        injection_count += n_states
        state_count += n_states
        parameter_count += n_params
        inner_vars_count += n_inner_vars
    end
    return wrapped_injector, state_count, parameter_count
end

function get_constant_lags(sys::PSY.System)
    delays = []
    for injector in get_injectors_with_dynamics(sys)
        device_delays = get_delays(PSY.get_dynamic_injector(injector))
        if !isnothing(device_delays)
            delays = vcat(delays, device_delays)
        end
    end
    return filter!(x -> x != 0, unique(delays))
    
end

function _wrap_dynamic_branches(
    sys::PSY.System,
    lookup::Dict{Int, Int},
    state_count::Int,
    parameter_count::Int,
)
    sys_base_power = PSY.get_base_power(sys)
    sys_base_freq = PSY.get_frequency(sys)
    dynamic_branches = get_dynamic_branches(sys)
    wrapped_branches = Vector{BranchWrapper}(undef, length(dynamic_branches))
    if !isempty(wrapped_branches)
        branches_count = 1
        for (ix, br) in enumerate(dynamic_branches)
            @debug "Wrapping Branch $(PSY.get_name(br))"
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = lookup[from_bus_number]
            bus_ix_to = lookup[to_bus_number]
            ix_range = range(state_count; length = n_states)
            ode_range = range(branches_count; length = n_states)
            n_params = get_n_params(br)
            p_range = range(parameter_count; length = n_params)
            @debug "ix_range=$ix_range ode_range=$ode_range p_range=$ode_range"
            wrapped_branches[ix] = BranchWrapper(
                br,
                bus_ix_from,
                bus_ix_to,
                ix_range,
                ode_range,
                p_range,
                sys_base_power,
                sys_base_freq,
            )
            branches_count += n_states
            state_count += n_states
            parameter_count += n_params
        end
    else
        @debug("System doesn't contain Dynamic Branches")
    end
    return wrapped_branches, state_count, parameter_count
end

function _wrap_static_injectors(
    sys::PSY.System,
    lookup::Dict{Int, Int},
    parameter_count::Int,
)
    static_injection_data = get_injection_without_dynamics(sys)
    container = Vector{StaticWrapper}(undef, length(static_injection_data))
    for (ix, ld) in enumerate(static_injection_data)
        if isa(ld, PSY.FixedAdmittance)
            continue
        end
        bus_n = PSY.get_number(PSY.get_bus(ld))
        bus_ix = lookup[bus_n]
        n_params = get_n_params(ld)
        p_range = range(parameter_count; length = n_params)
        container[ix] = StaticWrapper(ld, bus_ix, p_range)
        parameter_count += n_params
    end
    return container, parameter_count
end

function _wrap_loads(sys::PSY.System, lookup::Dict{Int, Int}, parameter_count::Int)
    sys_base_power = PSY.get_base_power(sys)
    # This needs to change if we implement dynamic load models
    static_loads =
        PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
    map_bus_load = Dict{PSY.Bus, Vector{PSY.ElectricLoad}}()
    for ld in static_loads
        if PSY.get_dynamic_injector(ld) !== nothing || !(PSY.get_available(ld))
            continue
        end
        bus = PSY.get_bus(ld)
        # Optimize this dictionary push
        push!(get!(map_bus_load, bus, PSY.ElectricLoad[]), ld)
    end
    return _construct_load_wrapper(lookup, map_bus_load, sys_base_power, parameter_count)
end

function _construct_load_wrapper(
    lookup::Dict{Int, Int},
    map_bus_load::Dict{PSY.Bus, Vector{PSY.ElectricLoad}},
    sys_base_power,
    parameter_count,
)
    container = Vector{StaticLoadWrapper}(undef, length(map_bus_load))
    for (ix, (bus, loads)) in enumerate(map_bus_load)
        bus_n = PSY.get_number(bus)
        bus_ix = lookup[bus_n]
        n_params = 8
        p_range = range(parameter_count; length = n_params)
        container[ix] = StaticLoadWrapper(bus, loads, p_range, bus_ix, sys_base_power)
        parameter_count += n_params
    end
    return container, parameter_count
end

function _get_ybus(sys::PSY.System)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    dyn_lines = PSY.get_components(x -> PSY.get_available(x), PSY.DynamicBranch, sys)
    if !isempty(PSY.get_components(PSY.ACBranch, sys))
        Ybus_ = PNM.Ybus(sys)
        ybus = Ybus_[:, :]
        lookup = Ybus_.lookup[1]
        ybus_rectangular = transform_ybus_to_rectangular(ybus)
        for br in dyn_lines
            ybus_update!(ybus_rectangular, br.branch, lookup, -1.0)
        end
    else
        ybus_rectangular =
            SparseArrays.SparseMatrixCSC{Float64, Int}(zeros(2 * n_buses, 2 * n_buses))
        Ybus_ = PNM.Ybus(sys)
        lookup = Ybus_.lookup[1]
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

function _get_shunt_values(br::PSY.Line)
    return PSY.get_b(br).from, PSY.get_b(br).to
end

function _get_shunt_values(::Union{PSY.TapTransformer, PSY.Transformer2W})
    return 0.0, 0.0
end

function _make_total_shunts(wrapped_branches, n_buses::Int)
    shunts = SparseArrays.SparseMatrixCSC{Float64, Int}(zeros(2 * n_buses, 2 * n_buses))
    for br in wrapped_branches
        bus_ix_from = get_bus_ix_from(br)
        bus_ix_to = get_bus_ix_to(br)
        b_from, b_to = _get_shunt_values(br.branch.branch)
        shunts[bus_ix_from, bus_ix_from + n_buses] += b_from
        shunts[bus_ix_to, bus_ix_to + n_buses] += b_to
        shunts[bus_ix_from + n_buses, bus_ix_from] -= b_from
        shunts[bus_ix_to + n_buses, bus_ix_to] -= b_to
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
    # Takes the upper quadrant of the rectangular shunts matrix
    shunts = LinearAlgebra.diag(total_shunts[1:n_buses, (n_buses + 1):end])
    for (ix, val) in enumerate(shunts)
        if val > 0
            @debug "Found shunt with value $val in bus index $ix"
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
) where {T <: Union{ConstantFrequency, ReferenceBus}}
    global_vars_dict = get_vars_ix()
    global_vars_dict[GLOBAL_VAR_SYS_FREQ_INDEX] = get_frequency_reference!(
        frequency_reference,
        wrapped_injectors,
        static_injection_data,
    )
    return global_vars_dict
end

function get_dynamic_wrapper(inputs::SimulationInputs, name::AbstractString)
    for w in get_dynamic_injectors(inputs)
        if PSY.get_name(w) == name
            return w
        end
    end
    println("DynamicWrapper with name $(name) was not found.")
    return nothing
end

function get_setpoints(inputs::SimulationInputs)
    dic = Dict{String, Dict{String, Float64}}()
    parameters = get_parameters(inputs)
    for w in get_dynamic_injectors(inputs)
        dic_w = Dict{String, Float64}()
        p_range = get_p_range(w)
        dic_w["P_ref"] = parameters[p_range][P_ref_ix]
        dic_w["Q_ref"] = parameters[p_range][Q_ref_ix]
        dic_w["ω_ref"] = parameters[p_range][ω_ref_ix]
        dic_w["V_ref"] = parameters[p_range][V_ref_ix]
        dic[PSY.get_name(w)] = dic_w
    end
    return dic
end
