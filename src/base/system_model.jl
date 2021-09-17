struct SystemModel{T <: SimulationModel, C <: Cache} <: Function
    inputs::SimulationInputs
    cache::C
end

function SystemModel{T}(inputs, cache::U) where {T <: SimulationModel, U <: Cache}
    return SystemModel{T, U}(inputs, cache)
end
