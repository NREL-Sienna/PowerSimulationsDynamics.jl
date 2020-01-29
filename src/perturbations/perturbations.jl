abstract type Perturbation end

struct ThreePhaseFault
    time::Float64
    Ybus #::SparseMatrixCSC{Complex{Float64}, Int64}
end

struct ControlStepUp
    time::Float64
    device::PSY.DynamicInjection
    signal::Symbol
    additional_value::Float64
end
