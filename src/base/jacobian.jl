# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

struct JacobianFunctionWrapper{
    F,
    T <: Union{Matrix{Float64}, SparseArrays.SparseMatrixCSC{Float64, Int64}},
} <: Function
    Jf::F
    Jv::T
    x::Vector{Float64}
end

function JacobianFunctionWrapper(
    m!::SystemModel{MassMatrixModel},
    x0_guess::Vector{Float64};
    # Improve the heuristic to do sparsity detection
    sparse_retrieve_loop = 3,
)
    x0 = deepcopy(x0_guess)
    n = length(x0)
    m_ = (residual, x) -> m!(residual, x, nothing, 0.0)
    jconfig = ForwardDiff.JacobianConfig(m_, similar(x0), x0, ForwardDiff.Chunk(x0))
    Jf = (Jv, x) -> begin
        @debug "Evaluating Jacobian Function"
        ForwardDiff.jacobian!(Jv, m_, zeros(n), x, jconfig)
        return
    end
    jac = zeros(n, n)
    for _ in 1:sparse_retrieve_loop
        temp = zeros(n, n)
        Jf(temp, abs.(x0 + Random.rand(n) - Random.rand(n)))
        jac .+= abs.(temp)
    end
    Jv = SparseArrays.sparse(jac)
    Jf(Jv, x0)
    return JacobianFunctionWrapper{typeof(Jf), typeof(Jv)}(Jf, Jv, x0)
end

function JacobianFunctionWrapper(
    m!::SystemModel{ResidualModel},
    x0::Vector{Float64};
    # Improve the heuristic to do sparsity detection
    sparse_retrieve_loop = 3,
)
    n = length(x0)
    m_ = (residual, x) -> m!(residual, zeros(n), x, nothing, 0.0)
    jconfig = ForwardDiff.JacobianConfig(m_, similar(x0), x0, ForwardDiff.Chunk(x0))
    Jf = (Jv, x) -> begin
        @debug "Evaluating Jacobian Function"
        ForwardDiff.jacobian!(Jv, m_, zeros(n), x, jconfig)
        return
    end
    jac = zeros(n, n)
    for _ in 1:sparse_retrieve_loop
        temp = zeros(n, n)
        Jf(temp, abs.(x0 + Random.rand(n) - Random.rand(n)))
        jac .+= abs.(temp)
    end
    Jv = SparseArrays.sparse(jac)
    Jf(Jv, x0)
    return JacobianFunctionWrapper{typeof(Jf), typeof(Jv)}(Jf, Jv, x0)
end

function (J::JacobianFunctionWrapper)(x::AbstractVector{Float64})
    J.x .= x
    return J.Jf(J.Jv, x)
end

function (J::JacobianFunctionWrapper)(
    JM::U,
    x::AbstractVector{Float64},
) where {U <: Union{Matrix{Float64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}}
    J.x .= x
    J.Jf(JM, x)
    return
end

function (J::JacobianFunctionWrapper)(
    JM::U,
    x::AbstractVector{Float64},
    p,
    t,
) where {U <: Union{Matrix{Float64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}}
    J.x .= x
    J.Jf(JM, x)
    return
end

function (J::JacobianFunctionWrapper)(
    JM::U,
    dx::AbstractVector{Float64},
    x::AbstractVector{Float64},
    p,
    gamma,
    t,
) where {U <: Union{Matrix{Float64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}}
    J.x .= x
    JM .= gamma * LinearAlgebra.Diagonal(ones(length(x))) .- J.Jf(JM, x)
    return JM
end

function get_jacobian(
    ::Type{T},
    inputs::SimulationInputs,
    x0_init::Vector{Float64},
) where {T <: SimulationModel}
    return JacobianFunctionWrapper(T(inputs, x0_init, JacobianCache), x0_init)
end

function _set_operating_point!(
    x0_init::Vector{Float64},
    inputs::SimulationInputs,
    system::PSY.System,
)
    status = power_flow_solution!(x0_init, system, inputs)
    status = initialize_static_injection!(inputs)
    status = initialize_dynamic_injection!(x0_init, inputs, system)
    status = initialize_dynamic_branches!(x0_init, inputs)
    return status
end
function get_jacobian(T, system::PSY.System)
    # Deepcopy avoid system modifications
    simulation_system = deepcopy(system)
    inputs = SimulationInputs(T, simulation_system, ReferenceBus)
    x0_init = get_flat_start(inputs)
    _set_operating_point!(x0_init, inputs, system)
    return get_jacobian(T, inputs, x0_init)
end
