"""
function f(::Simulation)

end 
g = get_parameter_sensitivity_function!(sim, Vector{Tuple{String, Type{T}, Symbol}}, f) where {T<: Union{Device, DynamicComponent} 
p = get_parameter_sensitivity_values(sim, Vector{Tuple{String, Type{T}, Symbol}}, f) where {T<: Union{Device, DynamicComponent} 

Zygote.gradient(g, p)
"""

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
    Symbol(s, :_OuterControl)
_append_symbol(s::Symbol, ::Type{T}) where {T <: PSY.InverterLimiter} =
    Symbol(s, :_InverterLimiter)

function get_indices_in_parameter_vector(sim, device_param_pairs)
    indices = Int[]
    for (device_name, component_type, param_symbol) in device_param_pairs
        ix = findfirst(
            x -> PSY.get_name(x.dynamic_device) == device_name,
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
                CRC.@ignore_derivatives @error "Device name not found in dynamic or static injectors"
            end
        end
        full_symbol = _append_symbol(param_symbol, component_type)
        external_ix = get_p_range(wrapped_device)
        internal_ix = findfirst(isequal(full_symbol), get_params_symbol(wrapped_device))
        global_ix = external_ix[internal_ix]
        @assert global_ix !== nothing
        push!(indices, global_ix)
        return indices
    end
end

function get_required_initialization_level(sim, device_param_pairs)
    return INITIALIZED #TODO -implement check of what type of parameters were changed
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

function get_parameter_sensitivity_function!(sim, device_param_pairs, f)
    indices = get_indices_in_parameter_vector(sim, device_param_pairs) 
    sim_level = get_required_initialization_level(sim, device_param_pairs)
    reset!(sim)
    @assert sim.status == BUILT
    sim.initialize_level = sim_level
    sim.enable_sensitivity = true
    sensitivity_function = (p) ->
        begin
            sim.inputs = deepcopy(sim.inputs_init)
            set_parameters!(sim, indices, p)   
            reset!(sim)   
            return f(sim)  
        end
    return sensitivity_function
end

function get_parameter_sensitivity_values(sim, device_param_pairs)
    indices = get_indices_in_parameter_vector(sim, device_param_pairs)
    param_vector = sim.inputs.parameters
    return param_vector[indices]
end
