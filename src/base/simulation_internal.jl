mutable struct SimulationInternal{T, U, V}
    x0_init::Vector{Float64}
    jacobian::T
    system_model::U
    sciml_function::V
    tstops::Vector{Float64}
    callbacks::DiffEqBase.CallbackSet
    problem::Union{Nothing, SciMLBase.DEProblem}
end

function SimulationInternal(x0_init, jacobian, system_model, model_function)
    return SimulationInternal(
        x0_init,
        jacobian,
        system_model,
        model_function,
        Vector{Float64}(),
        DiffEqBase.CallbackSet(),
        nothing,
    )
end

get_initial_conditions(internal::SimulationInternal) = internal.x0_init
get_jacobian(internal::SimulationInternal) = internal.jacobian
get_system_model(internal::SimulationInternal) = internal.system_model
