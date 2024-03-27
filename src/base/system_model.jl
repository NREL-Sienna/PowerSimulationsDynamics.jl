struct SystemModel{T <: SimulationModel, D <: DelayModel, C <: Cache}
    inputs::SimulationInputs
    cache::C
end

function SystemModel{T, D}(
    inputs,
    cache::U,
) where {T <: SimulationModel, D <: DelayModel, U <: Cache}
    return SystemModel{T, D, U}(inputs, cache)
end
