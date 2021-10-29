function update_global_vars!(
    cache::Cache,
    inputs::SimulationInputs,
    x::AbstractArray{U},
) where {U <: ACCEPTED_REAL_TYPES}
    for (var, pointer) in get_global_vars_update_pointers(inputs)
        pointer == 0 && continue
        get_global_vars(cache, U)[var] = x[pointer]
    end
    return
end
