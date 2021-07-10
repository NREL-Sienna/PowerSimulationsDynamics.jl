# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

struct JacobianFunctionWrapper{T}
    input::sim_inputs
    Jf::T
    Jv::SparseMatrixCSC{Float64, Int64}
    x::Vector{Float64}
    function modelJ(input, x)
        model = modelf(input, zeros(3))
        fu = (u) -> model(u, 0)
        jconfig = ForwardDiff.JacobianConfig(fu, x, ForwardDiff.Chunk{3}())
        Jf = (Jv, x) -> ForwardDiff.jacobian!(Jv, fu, x, jconfig)
        n = length(x)
        sJ = sparse(Jf(zeros(n, n), x))
        new{typeof(Jf)}(input, Jf, sJ, x)
    end
end

function (J::modelJ{T})(x) where {T <: Function}
    J.x .= x
    return J.Jf(J.Jv, x)
end
