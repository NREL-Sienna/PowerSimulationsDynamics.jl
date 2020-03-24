abstract type Perturbation end

struct ThreePhaseFault <: Perturbation
    time::Float64
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64}
end

get_affect(pert::ThreePhaseFault) =
    (integrator) -> PSY.get_ext(integrator.p[1])["Ybus"] = pert.Ybus

struct ControlReferenceChange <: Perturbation
    time::Float64
    device::PSY.DynamicInjection
    signal_index::Int64
    ref_value::Float64
end

function get_affect(pert::ControlReferenceChange)
    return (integrator) -> begin
        control_ref = PSY.get_ext(pert.device)[CONTROL_REFS]
        return control_ref[pert.signal_index] = pert.ref_value
    end
end
