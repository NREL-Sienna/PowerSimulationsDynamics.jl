"""
    function get_parameter_sensitivity_function!(
        sim::Simulation,
        device_parameter_pairs::Vector{Tuple{String, Type{T}, Symbol}},
        f::function,
    )

Gets a function for taking gradients with respect to parameters. 
# Arguments    
- `sim::Simulation` : Initialized simulation object
- `device_parameter_pairs::Vector{Tuple{String, Type{T}, Symbol}}` : Tuple used to identify the parameter, via the device name, as a `String`, the type of the Device or DynamicComponent, and the parameter as a `Symbol`. 
- `f::function` : User provided function with two inputs: a simulation and an additional input which can be used for data (```f(sim::Simulation, data::Any)```) The output must be a scalar value. This function can include executing the simulation and post-processing of results.  

# Example 
```julia
function f(sim::Simulation)
    execute!(sim, Rodas5())
    res = read_results(sim)
    _, δ = get_state_series(res, ("generator-1", :δ))
    sum(δ)
end 
g = get_parameter_sensitivity_function!(sim, ("generator-1", SingleMass, :H), f)
Zygote.gradient(g, [2.0]);
```
"""
function get_parameter_sensitivity_function!(sim, device_param_pairs, f)
    indices = get_indices_in_parameter_vector(sim, device_param_pairs)
    if indices === nothing
        return nothing
    end
    sim_level = get_required_initialization_level(sim, device_param_pairs)
    if sim_level === nothing
        return nothing
    end
    reset!(sim)
    @assert sim.status == BUILT
    sim.initialize_level = sim_level
    sim.enable_sensitivity = true
    sensitivity_function = (p, data) ->
        begin
            sim.inputs = deepcopy(sim.inputs_init)
            set_parameters!(sim, indices, p)
            reset!(sim)
            return f(sim, data)
        end
    return sensitivity_function
end

"""
    function get_parameter_sensitivity_function!(
        sim::Simulation,
        device_parameter_pairs::Vector{Tuple{String, Type{T}, Symbol}},
        f::function,
    )

get_parameter_sensitivity_values can be used in conjunction with get_parameter_sensitivity_function! to get the starting values of the parameters for taking gradients. 

# Arguments    
- `sim::Simulation` : Initialized simulation object
- `device_parameter_pairs::Vector{Tuple{String, Type{T}, Symbol}}` : Tuple used to identify the parameter, via the device name, as a `String`, the type of the Device or DynamicComponent, and the parameter as a `Symbol`. 

# Example 
```julia
function f(sim::Simulation)
    execute!(sim, Rodas5())
    res = read_results(sim)
    _, δ = get_state_series(res, ("generator-1", :δ))
    sum(δ)
end 
g = get_parameter_sensitivity_function!(sim, ("generator-1", SingleMass, :H), f)
p = get_parameter_sensitivity_values(sim, ("generator-1", SingleMass, :H))
Zygote.gradient(g, p);
```
"""
function get_parameter_sensitivity_values(sim, device_param_pairs)
    indices = get_indices_in_parameter_vector(sim, device_param_pairs)
    param_vector = sim.inputs.parameters
    return param_vector[indices]
end

function _append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.Device}
    return s
end
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.AVR} = Symbol(s, :_AVR)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.Machine} = Symbol(s, :_Machine)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.PSS} = Symbol(s, :_PSS)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.Shaft} = Symbol(s, :_Shaft)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.TurbineGov} = Symbol(s, :_TurbineGov)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.Converter} = Symbol(s, :_Converter)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.DCSource} = Symbol(s, :_DCSource)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.Filter} = Symbol(s, :_Filter)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.FrequencyEstimator} =
    Symbol(s, :_FrequencyEstimator)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.InnerControl} =
    Symbol(s, :_InnerControl)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.OuterControl} =
    @error "Specify PSY.ActivePowerControl or PSY.ReactivePowerControl"
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.ActivePowerControl} =
    Symbol(s, :_ActivePowerControl)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.ReactivePowerControl} =
    Symbol(s, :_ReactivePowerControl)

_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.InverterLimiter} =
    Symbol(s, :_InverterLimiter)

function get_indices_in_parameter_vector(sim, device_param_pairs)
    indices = Int[]
    for (device_name, component_type, param_symbol) in device_param_pairs
        ix = findfirst(
            x -> PSY.get_name(x.device) == device_name,
            sim.inputs_init.dynamic_injectors,
        )
        if ix !== nothing
            wrapped_device = sim.inputs_init.dynamic_injectors[ix]
        else
            ix = findfirst(
                x -> PSY.get_name(x.device) == device_name,
                sim.inputs_init.static_injectors,
            )
            if ix !== nothing
                wrapped_device = sim.inputs_init.static_injectors[ix]
            else
                CRC.@ignore_derivatives @warn "Device $device_name not found in dynamic or static injectors"
                return nothing
            end
        end
        full_symbol = _append_symbol(param_symbol, component_type)
        external_ix = get_p_range(wrapped_device)
        internal_ix = findfirst(isequal(full_symbol), get_params_symbol(wrapped_device))
        if internal_ix === nothing
            @warn "Parameter :$param_symbol of $component_type not found."
            return nothing
        end
        global_ix = external_ix[internal_ix]
        push!(indices, global_ix)
        return indices
    end
end

function get_required_initialization_level(sim, device_param_pairs)
    init_level = INITIALIZED
    for (device_name, component_type, param_symbol) in device_param_pairs
        metadata = get_params_metadata(component_type(nothing))
        symbols = [m.symbol for m in metadata]
        full_symbol = _append_symbol(param_symbol, component_type)
        ix = findfirst(isequal(full_symbol), symbols)
        metadata_entry = metadata[ix]
        if metadata_entry.in_mass_matrix == true
            @warn "Parameter :$param_symbol of $component_type appears in mass matrix -- not supported"
            return
        end
        if metadata_entry.in_network == true
            @warn "Parameter :$param_symbol of $component_type appears in network -- not supported"
            return
        end
        if metadata_entry.impacts_ic == true
            @warn "Parameter :$param_symbol of $component_type appears in initialization -- not supported"
            return
        end
        if metadata_entry.impacts_pf == true
            @warn "Parameter :$param_symbol of $component_type impacts power flow -- not supported"
            return
        end
    end
    return init_level
end

function make_buffer(a)
    buf = Zygote.Buffer(a)
    for i in eachindex(a)
        buf[i] = a[i]
    end
    return buf
end

function make_array(b)
    return copy(b)
end

#TODO - try to go back to making simulation inputs mutable and not using Accessors.jl?
#https://fluxml.ai/Zygote.jl/latest/limitations/#mutable-structs-1
#https://github.com/FluxML/Zygote.jl/issues/1127
function set_parameters!(sim, indices, params)
    inputs = sim.inputs
    parameter_buffer = make_buffer(inputs.parameters)
    for (ix, p) in zip(indices, params)
        parameter_buffer[ix] = p
    end
    Accessors.@reset inputs.parameters = make_array(parameter_buffer)
    sim.inputs = inputs
end
