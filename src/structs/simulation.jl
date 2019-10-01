mutable struct DynamicSimulation
    system::DynamicSystem
    problem::DiffEqBase.DAEProblem
    callbacks::Vector{DiffEqBase.DiscreteCallback}
    tstops::Vector{Float64}
    x0_init::Vector{Float64}
end

tstop = [1.0] #Define a timestop at t=1, the step change
cb = DiffEqBase.DiscreteCallback(LITS.change_t_one, LITS.Y_change!)

    prob = DiffEqBase.DAEProblem(system_model!, dx0, x0_init, tspan,
                            (Ybus_fault, case1_DynSystem), differential_vars = diff_vars)


function DynamicSimulation(system::DynamicSystem,
                            tspan,
