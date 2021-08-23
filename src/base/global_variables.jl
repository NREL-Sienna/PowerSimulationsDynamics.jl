function update_global_vars!(cache::Cache, inputs::SimulationInputs, x::AbstractArray{U}) where {U <: Real}
    index = get_global_vars_update_pointers(inputs)[GLOBAL_VAR_SYS_FREQ_INDEX]
    index == 0 && return
    # get_global_vars(cache)[:ω_sys] = x[index]
    return
end
