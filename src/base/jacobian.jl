# The implementation of caches is influenced by the SparseDiffTools.jl. We implement this
# custom code since the structure of the system models is not compatible with the functionalities
# in SparseDiffTools

struct JacobianFunctionWrapper{T}
    Jf::T
    Jv::SparseMatrixCSC{Float64, Int64}
    x::Vector{Float64}
    function JacobianFunctionWrapper(
        m!,
        x0;
        # Improve the heuristic to do sparsity detection
        sparse_retrieve_loop = 3,
    )
        n = length(x0)
        m_ = (residual, x) -> m!(residual, zeros(n), x, nothing, 0.0)
        jconfig = ForwardDiff.JacobianConfig(m_, similar(x0), x0, ForwardDiff.Chunk(x0))
        Jf = (Jv, x) -> ForwardDiff.jacobian!(Jv, fu, zeros(n), x, jconfig)
        jac_proto = zeros(n, n)
        for _ in 1:sparse_retrieve_loop
            temp = zeros(n, n)
            Jf(temp, abs(x0 + Random.rand(n) - Random.rand(n)))
            jac_proto .+= abs(temp)
        end
        Jv = sparse(jac_proto)
        new(Jf, Jv, x)
    end
end

function (J::JacobianFunctionWrapper{T})(x) where {T <: Function}
    J.x .= x
    return J.Jf(J.Jv, x)
end
