abstract type Perturbation end


mutable struct BranchTrip <: Perturbation
    time::Float64
    branch_name::String
end


function _ybus_update!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.ACBranch,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[PSY.get_number(arc.from)]
    bus_to_no = num_bus[arc.to.number]

    Y_l = (1 / (PSY.get_r(b) + PSY.get_x(b) * 1im))
    Y11 = Y_l + (1im * PSY.get_b(b).from)
    Y22 = Y_l + (1im * PSY.get_b(b).to)
    ybus[bus_from_no, bus_from_no] += mult * Y11
    Y12 = -Y_l
    ybus[bus_from_no, bus_to_no] += mult * Y12
    #Y21 = Y12
    ybus[bus_to_no, bus_from_no] += mult * Y12
    Y22 = Y_l + (1im * PSY.get_b(b).to)
    ybus[bus_to_no, bus_to_no] += mult * Y22
    return
end

function _ybus_update!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.TapTransformer,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[PSY.get_number(arc.from)]
    bus_to_no = num_bus[arc.to.number]

    Y_t = 1 / (PSY.get_r(b) + PSY.get_x(b) * 1im)
    c = 1 / PSY.get_tap(b)
    b = PSY.get_primary_shunt(b)

    Y11 = (Y_t * c^2)
    ybus[bus_from_no, bus_from_no] += mult * Y11
    Y12 = (-Y_t * c)
    ybus[bus_from_no, bus_to_no] += mult * Y12
    #Y21 = Y12
    ybus[bus_to_no, bus_from_no] += mult * Y12
    Y22 = Y_t
    ybus[bus_to_no, bus_to_no] += mult * (Y22 + (1im * b))

    return
end

function _ybus_update!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[PSY.get_number(arc.from)]
    bus_to_no = num_bus[arc.to.number]

    Y_t = 1 / (PSY.get_r(b) + PSY.get_x(b) * 1im)
    tap = (PSY.get_tap(b) * exp(PSY.get_α(b) * 1im))
    c_tap = (PSY.get_tap(b) * exp(-1 * PSY.get_α(b) * 1im))
    b = PSY.get_primary_shunt(b)

    Y11 = (Y_t / abs(tap)^2)
    ybus[bus_from_no, bus_from_no] += mult * Y11
    Y12 = (-Y_t / c_tap)
    ybus[bus_from_no, bus_to_no] += mult * Y12
    Y21 = (-Y_t / tap)
    ybus[bus_to_no, bus_from_no] += mult * Y21
    Y22 = Y_t
    ybus[bus_to_no, bus_to_no] += mult * (Y22 + (1im * b))
    return
end

function _ybus_update!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    fa::PSY.FixedAdmittance,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    bus = PSY.get_bus(fa)
    bus_ix = num_bus[PSY.get_number(bus)]
    ybus[bus_ix, bus_ix] += mult * fa.Y
    return
end


function get_affect(sys::PSY.System, pert::BranchTrip)
    branch = first(PSY.get_components(PSY.ACBranch,
                   sys, x -> PSY.get_name(x) == pert.branch_name)
                   )
    return (integrator) -> begin
        _ybus_update!(integrator.p.Ybus, branch, integrator.p.lookup, -1.0)
    end
    return

end

mutable struct BusTrip <: Perturbation
    time::Float64
    bus::PSY.Bus
end

function get_affect(::PSY.System, pert::BusTrip)
    return (integrator) -> begin
        0.0
    end
end

mutable struct ThreePhaseFault <: Perturbation
    time::Float64
    bus::PSY.Bus
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
