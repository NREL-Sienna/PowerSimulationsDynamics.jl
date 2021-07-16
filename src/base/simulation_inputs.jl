struct SimulationInputs
    base_power::Float64
    base_frequency::Float64
    injectors_data::Vector{DeviceWrapper{<:PSY.DynamicInjection}}
    static_injection_data::Vector{<:PSY.StaticInjection}
    dynamic_branches::Vector{BranchWrapper}
    injection_n_states::Int
    branches_n_states::Int
    variable_count::Int
    bus_count::Int
    ode_range::UnitRange{Int}
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}
    dyn_lines::Bool
    total_shunts::LinearAlgebra.Diagonal{Complex{Float64}}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::LinearAlgebra.Diagonal{Float64}
    tspan::NTuple{2, Float64}
    global_vars_update_pointers::Dict{Int, Vector{Int}}

    function SimulationInputs(sys::PSY.System, tspan::NTuple{2, Float64} = (0.0, 0.0))
        n_buses = get_n_buses(sys)
        Ybus, lookup = _get_Ybus(sys)
        wrapped_branches = _wrap_dynamic_branches(sys, lookup)
        has_dyn_lines = !isempty(wrapped_branches)
        branch_state_counts = 2 * length(wrapped_branches)
        injection_start = 2 * n_buses + branch_state_counts + 1
        wrapped_injectors = _wrap_injector_data(sys, lookup, injection_start)
        var_count = wrapped_injectors[end].ix_range[end]

        mass_matrix = _make_mass_matrix(wrapped_injectors, var_count, n_buses)
        DAE_vector = _make_DAE_vector(mass_matrix, var_count, n_buses)
        total_shunts = _make_total_shunts(wrapped_branches, n_buses)
        sys_f = PSY.get_frequency(sys)
        _adjust_states!(DAE_vector, mass_matrix, total_shunts, n_buses, sys_f)

        static_injection_data = PSY.get_components(
            PSY.StaticInjection,
            sys,
            x -> PSY.get_dynamic_injector(x) === nothing && PSY.get_available(x),
        )

        new(
            PSY.get_base_power(sys),
            sys_f,
            wrapped_injectors,
            collect(static_injection_data),
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
            tspan,
        )
    end
end

get_base_power(inputs::SimulationInputs) = inputs.base_power
get_base_frequency(inputs::SimulationInputs) = inputs.base_frequency
get_injectors_data(inputs::SimulationInputs) = inputs.injectors_data
get_dynamic_branches(inputs::SimulationInputs) = inputs.dynamic_branches
get_static_injections_data(inputs::SimulationInputs) = inputs.static_injection_data
get_voltage_buses_ix(inputs::SimulationInputs) = inputs.voltage_buses_ix
get_current_buses_ix(inputs::SimulationInputs) = inputs.current_buses_ix
get_Ybus(inputs::SimulationInputs) = inputs.Ybus
get_total_shunts(inputs::SimulationInputs) = inputs.total_shunts
get_lookup(inputs::SimulationInputs) = inputs.lookup
get_DAE_vector(inputs::SimulationInputs) = inputs.DAE_vector
get_mass_matrix(inputs::SimulationInputs) = inputs.mass_matrix
get_tspan(inputs::SimulationInputs) = inputs.tspan
has_dyn_lines(inputs::SimulationInputs) = inputs.dyn_lines
get_global_vars_update_pointers(inputs::SimulationInputs) =
    inputs.global_vars_update_pointers

get_injection_n_states(inputs::SimulationInputs) = inputs.injection_n_states
get_branches_n_states(inputs::SimulationInputs) = inputs.branches_n_states
get_variable_count(inputs::SimulationInputs) = inputs.variable_count
get_ode_range(inputs::SimulationInputs) = inputs.ode_range
get_bus_count(inputs::SimulationInputs) = inputs.bus_count
get_bus_range(inputs::SimulationInputs) = 1:(2 * inputs.bus_count)
# get_ω_sys(inputs::SimulationInputs) = get_global_vars(inputs)[:ω_sys]

function _wrap_injector_data(sys::PSY.System, lookup, injection_start::Int)
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
            DeviceWrapper(device, bus_ix, ix_range, ode_range, inner_vars_range)
        injection_count += n_states
        injection_start += n_states
        inner_vars_count = inner_vars_range[end]
    end
    return wrapped_injector
end

function _wrap_dynamic_branches(sys::PSY.System, lookup::Dict{Int, Int})
    branches_start = 2*get_n_buses(sys) + 1
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
    total_shunts = LinearAlgebra.Diagonal(zeros(Complex{Float64}, n_buses))
    if isempty(wrapped_branches)
        return total_shunts
    else
        for br in wrapped_branches
            bus_ix_from = get_bus_ix_from(br)
            bus_ix_to = get_bus_ix_to(br)
            b_from = PSY.get_b(br).from
            b_to = PSY.get_b(br).to
            total_shunts[bus_ix_from, bus_ix_from] += 1im * b_from
            total_shunts[bus_ix_to, bus_ix_to] += 1im * b_to
        end
    end
    return total_shunts
end

function _adjust_states!(
    DAE_vector::BitVector,
    mass_matrix::LinearAlgebra.Diagonal{Float64},
    total_shunts::LinearAlgebra.Diagonal{Complex{Float64}},
    n_buses::Int,
    sys_f::Float64,
)
    all(iszero.(total_shunts)) && return
    line_constant = 1 / (2.0 * π * sys_f)
    for (ix, shunt) in enumerate(total_shunts)
        val = imag(shunt)
        if val > 0
            mass_matrix[ix, ix] =
                mass_matrix[ix + n_buses, ix + n_buses] = val * line_constant
            DAE_vector[ix] = DAE_vector[ix + n_buses] = true
        end
    end
    return
end

# Default implementation for both models. This implementation is to future proof if there is
# a divergence between the required build methods
function _build!(inputs::SimulationInputs, sys::PSY.System)
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
