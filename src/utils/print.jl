function Base.show(io::IO, op_model::PSY.System)
    println(io, "System()")
end

function Base.show(io::IO, op_model::DynamicSimulation)
    println(io, "Simulation()")
end
