#TODO: Make inmmutable later. Requires refactor of ThreePhase Fault callbacks
mutable struct SimulationInputs <: SimulationInputs
    sys::PSY.System
    injectors_data::Vector{<:PSY.StaticInjection}
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
    dynamic_injectors = Vector{PSY.StaticInjection}(),
    counts::Base.ImmutableDict{Symbol, Int} = Base.ImmutableDict{Symbol, Int}(),
    Ybus::SparseMatrixCSC{Complex{Float64}, Int} = SparseMatrixCSC{
        Complex{Float64},
        Int,
    }(
        zeros(1, 1),
    ),
    dyn_lines::Bool = false,
    voltage_buses_ix::Vector{Int} = Vector{Int}(),
    current_buses_ix::Vector{Int} = Vector{Int}(),
    global_index::Dict{String, Dict{Symbol, Int}} = Dict{String, Dict{Symbol, Int}}(),
    total_shunts::Dict{Int, Float64} = Dict{Int, Float64}(),
    global_vars::Dict{Symbol, Number} = Dict{Symbol, Number}(),
    lookup::Dict{Int, Int} = Dict{Int, Int}(),
    DAE_vector::Vector{Bool} = Vector{Bool}(),
    aux_arrays::Dict{Int, Vector} = Dict{Int, Vector}(),
    tspan::NTuple{2, Float64} = (0.0, 0.0),
)
    return SimulationInputs(
        sys,
        dynamic_injectors,
        counts,
        Ybus,
        dyn_lines,
        voltage_buses_ix,
        current_buses_ix,
        global_index,
        total_shunts,
        global_vars,
        lookup,
        DAE_vector,
        aux_arrays,
        tspan,
    )
end

function _add_dynamic_bus_states!(
    DAE_vector::Vector{Bool},
    voltage_buses_ix::Vector{Int},
    bus_ix::Int,
    n_buses::Int,
)
    push!(voltage_buses_ix, bus_ix)
    DAE_vector[bus_ix] = DAE_vector[bus_ix + n_buses] = true
    return
end

function _add_to_total_shunts!(total_shunts::Dict{Int, Float64}, pairs...)
    merge!(+, total_shunts, Dict(pairs...))
    return
end

function _index_dynamic_lines!(
    inputs::SimulationInputs,
    branch::PSY.DynamicBranch,
    n_buses::Int,
)
    DAE_vector = get_DAE_vector(inputs)
    voltage_buses_ix = get_voltage_buses_ix(inputs)
    arc = PSY.get_arc(branch)
    from_bus_number = PSY.get_number(arc.from)
    to_bus_number = PSY.get_number(arc.to)
    bus_ix_from = get_lookup(inputs)[from_bus_number]
    bus_ix_to = get_lookup(inputs)[to_bus_number]
    b_from = PSY.get_b(branch).from
    b_to = PSY.get_b(branch).to
    total_shunts = get_total_shunts(inputs)
    b_from > 0.0 && _add_to_total_shunts!(total_shunts, bus_ix_from => b_from)
    b_to > 0.0 && _add_to_total_shunts!(total_shunts, bus_ix_to => b_to)
    b_from > 0.0 &&
        _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_from, n_buses)
    b_to > 0.0 && _add_dynamic_bus_states!(DAE_vector, voltage_buses_ix, bus_ix_to, n_buses)
    n_states = PSY.get_n_states(branch)
    DAE_vector = push!(DAE_vector, collect(trues(n_states))...)
    return
end

function build!(inputs::SimulationInputs)
    sys = get_system(inputs)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    DAE_vector = inputs.DAE_vector = collect(falses(n_buses * 2))
    global_state_index = inputs.global_index = MAPPING_DICT()
    state_space_ix = [n_buses * 2]
    current_buses_ix = inputs.current_buses_ix = collect(1:n_buses)
    static_bus_var_count = 2 * n_buses
    voltage_buses_ix = inputs.voltage_buses_ix = Vector{Int}()
    total_states = 0
    first_dyn_branch_point = -1
    branches_n_states = 0
    global_vars =
        inputs.global_vars = Dict{Symbol, Number}(
            :ω_sys => 1.0,
            :ω_sys_index => -1, #To define 0 if infinite source, bus_number otherwise,
        )
    inputs.total_shunts = Dict{Int, Float64}()
    found_ref_bus = false

    #jd/TODO: Check logic here
    inputs.Ybus, inputs.lookup = _get_Ybus(sys)
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sys)

    if !(isempty(dyn_branches))
        inputs.dyn_lines = true
        first_dyn_branch_point = state_space_ix[1] + 1
        for br in dyn_branches
            _index_dynamic_lines!(inputs, br, n_buses)
            total_states += PSY.get_n_states(br)
            _add_states_to_global!(global_state_index, state_space_ix, br)
        end

        for (ix, val) in enumerate(DAE_vector[1:n_buses])
            if val
                #This V_ix should be V_number.
                global_state_index["V_$(ix)"] = Dict(:R => ix, :I => ix + n_buses)
                total_states += 2
                static_bus_var_count -= 2
                push!(voltage_buses_ix, ix)
                @assert static_bus_var_count >= 0
            end
        end
        branches_n_states = state_space_ix[1] - n_buses * 2
    else
        inputs.dyn_lines = false
        @debug("System doesn't contain Dynamic Branches")
    end
    unique!(voltage_buses_ix)
    sources = PSY.get_components(PSY.Source, sys)
    for s in sources
        btype = PSY.get_bustype(PSY.get_bus(s))
        if (btype == PSY.BusTypes.REF) && found_ref_bus
            throw(
                IS.ConflictingInputsError(
                    "The system can't have more than one source or generator in the REF Bus",
                ),
            )
        end
        btype != PSY.BusTypes.REF && continue
        global_vars[:ω_sys_index] = 0 #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end

    dynamic_injection = PSY.get_components(
        PSY.StaticInjection,
        sys,
        x -> PSY.get_dynamic_injector(x) !== nothing,
    )

    if isempty(dynamic_injection)
        error("System doesn't contain any DynamicInjection devices")
    end

    for d in dynamic_injection
        @debug PSY.get_name(d)
        dynamic_device = PSY.get_dynamic_injector(d)
        isempty(PSY.get_states(dynamic_device)) && continue
        device_bus = PSY.get_bus(d)
        btype = PSY.get_bustype(device_bus)
        if (btype == PSY.BusTypes.REF) && found_ref_bus
            throw(
                IS.ConflictingInputsError(
                    "The system can't have more than one source or generator in the REF Bus",
                ),
            )
        end
        state_types = make_device_index!(d)
        device_n_states = PSY.get_n_states(dynamic_device)
        DAE_vector = push!(DAE_vector, state_types...)
        total_states += device_n_states
        _add_states_to_global!(global_state_index, state_space_ix, dynamic_device)
        push!(inputs.injectors_data, d)

        btype != PSY.BusTypes.REF && continue
        if typeof(dynamic_device) <: PSY.DynamicGenerator
            ω_ix = global_state_index[PSY.get_name(d)][:ω]
        elseif typeof(dynamic_device) <: PSY.DynamicInverter
            #TO DO: Make it general for cases when ω is not a state (droop)!
            ω_ix = global_state_index[PSY.get_name(d)][:ω_oc]
        end
        global_vars[:ω_sys_index] = ω_ix #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end
    injection_n_states = state_space_ix[1] - branches_n_states - n_buses * 2
    @assert total_states == state_space_ix[1] - static_bus_var_count
    @debug total_states
    setdiff!(current_buses_ix, voltage_buses_ix)

    inputs.counts = Base.ImmutableDict(
        :total_states => total_states,
        :injection_n_states => injection_n_states,
        :branches_n_states => branches_n_states,
        :first_dyn_injection_pointer => 2 * n_buses + branches_n_states + 1,
        :first_dyn_branch_point => first_dyn_branch_point,
        :total_variables => total_states + static_bus_var_count,
        :bus_count => n_buses,
    )

    @assert get_ω_sys(inputs) != -1
    return
end

get_system(inputs::SimulationInputs) = inputs.sys
get_injectors_data(inputs::SimulationInputs) = inputs.injectors_data
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
