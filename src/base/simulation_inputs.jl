struct SimulationInputs
    base_power::Float64
    base_frequency::Float64
    injectors_data::Vector{DeviceWrapper{<:PSY.DynamicInjection}}
    static_injection_data::Vector{<:PSY.StaticInjection}
    dynamic_branches::Vector{BranchWrapper}
    injection_n_states::Int
    branches_n_states::Int
    total_variables::Int #var_count,
    bus_count::Int # n_buses,
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}
    dyn_lines::Bool
    voltage_buses_ix::Vector{Int}
    current_buses_ix::Vector{Int}
    total_shunts::LinearAlgebra.Diagonal{Complex{Float64}}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    mass_matrix::LinearAlgebra.Diagonal{Float64}
    tspan::NTuple{2, Float64}
    global_vars_update_pointers::Dict{Int, Vector{Int}}

    function SimulationInputs(sys::PSY.System, tspan::NTuple{2, Float64} = (0.0, 0.0))
        n_buses = length(PSY.get_bus_numbers(sys))
        Ybus, lookup = _get_Ybus(sys)
        dynamic_branches =
            collect(PSY.get_components(PSY.DynamicBranch, sys, x -> PSY.get_available(x)))
        branch_state_counts = 2 * length(dynamic_branches)

        wrapped_injectors = _wrap_injector_data(sys, lookup, 2 * n_buses + branch_state_counts + 1)


        var_count = wrapped_injectors[end].ix_range[end]

        mass_matrix = LinearAlgebra.Diagonal(ones(10))
        mass_matrix[1:(2 * n_buses), 1:(2 * n_buses)] .= 0.0

        static_injection_data = PSY.get_components(
            PSY.StaticInjection,
            sys,
            x -> PSY.get_dynamic_injector(x) === nothing && PSY.get_available(x),
        )

        new(
            PSY.get_base_power(sys),
            PSY.get_frequency(sys),
            wrapped_injectors,
            collect(static_injection_data),
            Vector{BranchWrapper}(),
            0,
            0,
            0,
            0,
            Ybus,
            !isempty(dynamic_branches),
            Vector{Int}(),
            collect(1:n_buses),
            LinearAlgebra.Diagonal(zeros(Complex{Float64}, n_buses)),
            lookup,
            _init_DAE_vector!(var_count, n_buses),
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
get_global_vars_update_pointers(inputs::SimulationInputs) = inputs.global_vars_update_pointers

# get_injection_pointer(inputs::SimulationInputs) =get_counts(inputs)[:first_dyn_injection_pointer]
# get_branches_pointer(inputs::SimulationInputs) = get_counts(inputs)[:first_dyn_branch_point]
get_n_injection_states(inputs::SimulationInputs) = inputs.injection_n_states
get_n_branches_states(inputs::SimulationInputs) = inputs.branches_n_states
get_variable_count(inputs::SimulationInputs) = inputs.total_variables
get_device_index(inputs::SimulationInputs, device::D) where {D <: PSY.DynamicInjection} =
    get_global_index(inputs)[device.name]
get_bus_count(inputs::SimulationInputs) = inputs.bus_count
# get_ω_sys(inputs::SimulationInputs) = get_global_vars(inputs)[:ω_sys]


function _wrap_injector_data(sys, lookup, injection_start)

    injector_data = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing && PSY.get_available(x),
    )

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
        wrapped_injector[ix] = DeviceWrapper(device, bus_ix, ix_range, ode_range, inner_vars_range)
        injection_count += n_states
        injection_start += n_states
        inner_vars_count = inner_vars_range[end]
    end
    return wrapped_injector
end

function _wrap_dynamic_branches(sys, lookup)
    branches_count = 1
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

function make_global_state(wrapped_devices
)
    global_state_index[PSY.get_name(device)] = Dict{Symbol, Int}()
    for s in PSY.get_states(device)
        state_space_ix[1] += 1
        global_state_index[PSY.get_name(device)][s] = state_space_ix[1]
    end

    return
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
        # Temporary condition needs to be removed for other static injections
        PSY.get_dynamic_injector(s) !== nothing && continue
        index_static_injection(inputs, s)
        _attach_control_refs!(s)
    end
    return
end


function _mass_matrix_inputs!(inputs::SimulationInputs)
    for d in PSID.get_injectors_data(inputs)
        device_mass_matrix_entries!(inputs, d.device)
    end
    if has_dyn_lines(inputs)
        mass_matrix = get_mass_matrix(inputs)
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

function _dae_vector_update!(inputs::SimulationInputs)
    mass_matrix = get_mass_matrix(inputs)
    DAE_vector = get_DAE_vector(inputs)
    bus_vars_count = 2 * get_bus_count(inputs)
    bus_range = 1:bus_vars_count
    injection_start = get_injection_pointer(inputs)
    ode_range = injection_start:get_variable_count(inputs)
    mass_buses = @view mass_matrix[bus_range, bus_range]
    DAE_ode = @view DAE_vector[ode_range]
    mass_ode = @view mass_matrix[ode_range, ode_range]
    for i in eachindex(DAE_vector[bus_range])
        IS.@assert_op DAE_vector[bus_range][i] == (mass_buses[i, i] > 0.0)
    end
    for i in eachindex(DAE_ode)
        DAE_ode[i] = (mass_ode[i, i] > 0.0)
    end
    if has_dyn_lines(inputs)
        branches_range =
            range(get_branches_pointer(inputs), length = get_n_branches_states(inputs))
        mass_branches = @view mass_matrix[branches_range, branches_range]
        for i in eachindex(DAE_vector[branches_range])
            IS.@assert_op DAE_vector[branches_range][i] == (mass_branches[i, i] > 0.0)
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

    #set_frequency_reference!(inputs, sys)
    #IS.@assert_op get_ω_sys(inputs) != -1

    _mass_matrix_inputs!(inputs)
    _dae_vector_update!(inputs)

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
