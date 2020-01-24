function Base.show(io::IO, op_model::PSY.System)
    println(io, "System()")
end

function Base.show(io::IO, op_model::Simulation)
    println(io, "Simulation()")
end
