function Base.show(io::IO, ::Simulation)
    println(io, "Simulation()")
end

function Base.show(io::IO, ::SimulationInputs)
    println(io, "SimulationInputs()")
end

function Base.show(io::IO, smr::SmallSignalOutput)
    val = smr.stable ? "is" : "is not"
    println(io, "The system $(val) small signal stable")
end

function Base.show(io::IO, pert::ControlReferenceChange)
    println(
        io,
        "Change of $(pert.signal) at time t = $(pert.time) of device $(PSY.get_name(pert.device)) to $(pert.ref_value) per unit",
    )
end

function Base.show(io::IO, ::MIME"text/html", res::SimulationResults)
    show_results_table(io, res, backend = :html)
    println(io)
end

function Base.show(io::IO, ::MIME"text/plain", res::SimulationResults)
    show_results_table(io, res, backend = :auto)
    println(io)
end

function show_results_table(io::IO, res::SimulationResults; kwargs...)
    header = ["Property", "Value"]
    sys = get_system(res)
    sol = get_solution(res)
    table = [
        "System Base Power [MVA]" string(PSY.get_base_power(sys))
        "System Base Frequency [Hz]" string(PSY.get_frequency(sys))
        "Time Span" string((sol.t[1], sol.t[end]))
        "Total Time Steps" string(length(sol.t))
        "Number of States" string(length(sol.u[1]))
    ]
    PrettyTables.pretty_table(
        io,
        table,
        header;
        title = "Simulation Results Summary",
        alignment = :l,
        kwargs...,
    )
end
