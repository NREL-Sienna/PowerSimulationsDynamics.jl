struct SmallSignalOutput
    jacobian::Matrix{Float64}
    eigenvalues::Vector{Complex{Float64}}
    eigenvectors::Matrix{Complex{Float64}}
    stable::Bool
    operating_point::Vector{Float64}
end
