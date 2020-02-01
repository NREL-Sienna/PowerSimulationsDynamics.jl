function Base.show(io::IO, ::Simulation)
    println(io, "Simulation()")
end

function Base.show(io::IO, smr::SmallSignalOutput)
    val = smr.stable ? "is" : "is not"
    println(io, "The system $(val) small signal stable")
end
