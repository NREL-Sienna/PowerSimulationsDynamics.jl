function Base.show(io::IO, op_model::DynamicSystem)
    println(io, "System()")
end

function Base.show(io::IO, op_model::DynamicSimulation)
    println(io, "Simulation()")
end
