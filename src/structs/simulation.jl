mutable struct DynamicSimulation
    dyn_system::DynamicSystem
    problem::DiffEqBase.DAEProblem
    callbacks::DiffEqBase.DiscreteCallback
    tstops::Vector{Float64}
    x0_init::Vector{Float64}
    initialized::Bool
    solution::Union{Nothing, DiffEqBase.DAESolution}
end

function DynamicSimulation(dyn_system::DynamicSystem,
                           tspan,
                           control,
                           callback,
                           x0_init)

    if !is_indexed(dyn_system)
        _index_dynamic_system!(dyn_system)
    end

    dx0 = zeros(get_total_rows(dyn_system))
    prob = DiffEqBase.DAEProblem(system_model!,
                              dx0,
                              x0_init,
                              tspan,
                              (control, dyn_system),
                              differential_vars = dyn_system.DAE_vector)

    return DynamicSimulation(dyn_system,
                             prob,
                             callback,
                             [1.0],
                             x0_init,
                             true,
                             nothing)

end


function run_simulation!(sim::DynamicSimulation, solver; kwargs...)
    sim.solution = DiffEqBase.solve(sim.problem,
                           solver;
                           callback = sim.callbacks,
                           tstops = sim.tstops, kwargs...)

    return
end
