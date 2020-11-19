abstract type Perturbation end


mutable struct BranchTrip <: Perturbation
    time::Float64
    branch::ACBranch
end

function get_affect(::PSY.System, pert::BranchTrip)
    return (integrator) -> begin
        0.0
    end
end

mutable struct BusTrip <: Perturbation
    time::Float64
    bus::Bus
end

function get_affect(::PSY.System, pert::BranchTrip)
    return (integrator) -> begin
        0.0
    end
end

mutable struct ThreePhaseFault <: Perturbation
    time::Float64
    bus::Bus
    impedance::Complex
end

function get_affect(::PSY.System, pert::ThreePhaseFault)
    return (integrator) -> begin
        0.0
    end
end

mutable struct NetworkSwitch <: Perturbation
    time::Float64
    Ybus::SparseMatrixCSC{Complex{Float64}, Int64}
end

function get_affect(::PSY.System, pert::NetworkSwitch)
    return (integrator) -> begin
        integrator.p.Ybus = pert.Ybus
    end
end

mutable struct ControlReferenceChange <: Perturbation
    time::Float64
    device::PSY.DynamicInjection
    signal_index::Int64
    ref_value::Float64
end

function get_affect(system::PSY.System, pert::ControlReferenceChange)
    device = PSY.get_component(typeof(pert.device), system, PSY.get_name(pert.device))
    pert.device = device
    return (integrator) -> begin
        control_ref = PSY.get_ext(device)[CONTROL_REFS]
        return control_ref[pert.signal_index] = pert.ref_value
    end
end
