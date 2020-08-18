#TODO: Make inmmutable later. Requires refactor of ThreePhase Fault callbacks
mutable struct SimulationInputs
    sys::PSY.System
    counts::Base.ImmutableDict{Symbol, Int}
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64}
    dyn_lines::Bool
    voltage_buses_ix::Vector{Int}
    current_buses_ix::Vector{Int}
    global_index::Dict{String, Dict{Symbol, Int64}}
    total_shunts::Dict{Int64, Float64}
    global_vars::Dict{Symbol, Number}
    lookup::Dict{Int, Int}
    DAE_vector::Vector{Bool}
    aux_arrays::Dict{Int, Vector}
    tspan::NTuple{2, Float64}
end

function SimulationInputs(;
    sys::PSY.System,
    counts::Base.ImmutableDict{Symbol, Int} = Base.ImmutableDict{Symbol, Int}(),
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64} = SparseMatrixCSC{Complex{Float64}, Int64}(zeros(1, 1)),
    dyn_lines::Bool = false,
    voltage_buses_ix::Vector{Int} = Vector{Int}(),
    current_buses_ix::Vector{Int} = Vector{Int}(),
    global_index::Dict{String, Dict{Symbol, Int64}} = Dict{String, Dict{Symbol, Int64}}(),
    total_shunts::Dict{Int64, Float64} = Dict{Int64, Float64}(),
    global_vars::Dict{Symbol, Number} = Dict{Symbol, Number}(),
    lookup::Dict{Int, Int} = Dict{Int, Int}(),
    DAE_vector::Vector{Bool} = Vector{Bool}(),
    aux_arrays::Dict{Int, Vector} = Dict{Int, Vector}(),
    tspan::NTuple{2, Float64} = (0.0, 0.0),
)
    return SimulationInputs(
        sys,
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

function build(inputs::SimulationInputs)
    sys = get_system(inputs)
    n_buses = length(PSY.get_components(PSY.Bus, sys))
    DAE_vector = collect(falses(n_buses * 2))
    global_state_index = MAPPING_DICT()
    state_space_ix = [n_buses * 2]
    current_buses_ix = collect(1:n_buses)
    static_bus_var_count = 2 * length(current_buses_ix)
    voltage_buses_ix = Vector{Int}()
    total_states = 0
    first_dyn_branch_point = -1
    branches_n_states = 0
    global_vars = Dict{Symbol, Number}(
        :ω_sys => 1.0,
        :ω_sys_index => -1, #To define 0 if infinite source, bus_number otherwise,
    )
    total_shunts = Dict{Int, Float64}()
    found_ref_bus = false
    sys_basepower = PSY.get_base_power(sys)

    inputs.Ybus, inputs.lookup = _get_Ybus(sys)
    dyn_branches = PSY.get_components(PSY.DynamicBranch, sys)

    if !(isempty(dyn_branches))
        inputs.dyn_lines = true
        first_dyn_branch_point = state_space_ix[1] + 1
        for br in dyn_branches
            arc = PSY.get_arc(br)
            n_states = PSY.get_n_states(br)
            from_bus_number = PSY.get_number(arc.from)
            to_bus_number = PSY.get_number(arc.to)
            bus_ix_from = lookup[from_bus_number]
            bus_ix_to = lookup[to_bus_number]
            merge!(
                +,
                total_shunts,
                Dict(
                    bus_ix_from => 1 / PSY.get_b(br).from,
                    bus_ix_to => 1 / PSY.get_b(br).to,
                ),
            )
            push!(voltage_buses_ix, bus_ix_from, bus_ix_to)
            DAE_vector[bus_ix_from] = DAE_vector[bus_ix_from + n_buses] = true
            DAE_vector[bus_ix_to] = DAE_vector[bus_ix_to + n_buses] = true
            DAE_vector = push!(DAE_vector, collect(trues(n_states))...)
            total_states += n_states
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
            throw(IS.ConflictingInputsError("The system can't have more than one source or generator in the REF Bus"))
        end
        btype != PSY.BusTypes.REF && continue
        global_vars[:ω_sys_index] = 0 #To define 0 if infinite source, bus_number otherwise,
        found_ref_bus = true
    end
    dynamic_injection = PSY.get_components(PSY.DynamicInjection, sys)
    isempty(dynamic_injection) &&
        error("System doesn't contain any DynamicInjection devices")
    for d in dynamic_injection
        if !(:states in fieldnames(typeof(d)))
            continue
        end
        device_bus = PSY.get_bus(d)
        btype = PSY.get_bustype(device_bus)
        if (btype == PSY.BusTypes.REF) && found_ref_bus
            throw(IS.ConflictingInputsError("The system can't have more than one source or generator in the REF Bus"))
        end
        _make_device_index!(d)
        device_n_states = PSY.get_n_states(d)
        DAE_vector = push!(DAE_vector, collect(trues(device_n_states))...)
        total_states += device_n_states
        _add_states_to_global!(global_state_index, state_space_ix, d)
        btype != PSY.BusTypes.REF && continue
        if typeof(d) <: PSY.DynamicGenerator
            ω_ix = global_state_index[PSY.get_name(d)][:ω]
        elseif typeof(d) <: PSY.DynamicInverter
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

    inputs.DAE_vector = DAE_vector
    inputs.global_index = global_state_index
    inputs.voltage_buses_ix = voltage_buses_ix
    inputs.current_buses_ix = current_buses_ix
    inputs.total_shunts = total_shunts
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
