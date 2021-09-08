abstract type Perturbation end

"""
Use to model the trip of an AC Branch in in the system. Accepts any ACBranch
"""
mutable struct BranchTrip <: Perturbation
    time::Float64
    branch_type::Type{<:PSY.ACBranch}
    branch_name::String
end

"""
Use to model changes in the network as sudden impedance changes. Mostly used for small system analysis.
"""
mutable struct BranchImpedanceChange <: Perturbation
    time::Float64
    branch_type::Type{<:PSY.ACBranch}
    branch_name::String
    multiplier::Float64
end

function _get_branch_for_perturbation(
    sys::PSY.System,
    perturbation::T,
) where {T <: Perturbation}
    if perturbation.branch_type == PSY.DynamicBranch
        error("DynamicBranch is not supported currently with perturbation $T")
    end
    branch = PSY.get_component(perturbation.branch_type, sys, perturbation.branch_name)
    if branch === nothing || !PSY.get_available(branch)
        throw(
            IS.ConflictingInputsError(
                "Branch $(perturbation.branch_type) $(perturbation.branch_name) not existent or unavailable in the system. Use PSY.show_components(system, $(perturbation.branch_type)) to check valid component names",
            ),
        )
    end
    return branch
end

function get_affect(sys::PSY.System, pert::BranchImpedanceChange)
    branch = _get_branch_for_perturbation(sys, pert)
    mult = 0.0
    if pert.multiplier < 0.0
        throw(
            IS.ConflictingInputsError(
                "Negative multipliers are not allowed for BranchImpedanceChange perturbations",
            ),
        )
    elseif pert.multiplier < 1.0
        mult = 1.0 - pert.multiplier
    elseif pert.multiplier > 1.0
        mult = pert.multiplier - 1.0
    else
        @assert pert.multiplier == 1.0
    end

    return (integrator) -> begin
        ybus_update!(integrator.p, branch, mult)
    end

    return
end

function get_affect(sys::PSY.System, pert::BranchTrip)
    branch = _get_branch_for_perturbation(sys, pert)
    return (integrator) -> begin
        ybus_update!(integrator.p, branch, -1.0)
    end
    return
end

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    b::PSY.ACBranch,
    num_bus::Dict{Int, Int},
    mult::Float64,
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

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    b::PSY.TapTransformer,
    num_bus::Dict{Int, Int},
    mult::Float64,
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

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    b::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int, Int},
    mult::Float64,
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

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{ComplexF64, Int},
    fa::PSY.FixedAdmittance,
    num_bus::Dict{Int, Int},
    mult::Float64,
)
    bus = PSY.get_bus(fa)
    bus_ix = num_bus[PSY.get_number(bus)]
    ybus[bus_ix, bus_ix] += mult * fa.Y
    return
end

function ybus_update!(integrator_params, branch::PSY.ACBranch, mult::Float64)
    ybus_update!(integrator_params.ybus_rectangular, branch, integrator_params.lookup, mult)
    return
end

# Pending implementation
#=
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
=#

"""
Use to model a change in the network by switching the underlying Ybus in the simulation
"""
mutable struct NetworkSwitch <: Perturbation
    time::Float64
    Ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int}
end

function get_affect(::PSY.System, pert::NetworkSwitch)
    return (integrator) -> begin
        # TODO: This code can be more performant using SparseMatrix methods
        for (i, v) in enumerate(pert.Ybus)
            integrator.p.Ybus[i] = v
        end
    end
end

"""
Use to model control reference changes in devices of the model
"""
mutable struct ControlReferenceChange <: Perturbation
    time::Float64
    device::PSY.DynamicInjection
    signal_index::Int
    ref_value::Float64
end

# TODO: change this to use setters
function get_affect(system::PSY.System, pert::ControlReferenceChange)
    device = PSY.get_component(typeof(pert.device), system, PSY.get_name(pert.device))
    pert.device = device
    return (integrator) -> begin
        control_ref = PSY.get_ext(device)[CONTROL_REFS]
        return control_ref[pert.signal_index] = pert.ref_value
    end
end

"""
Use to model a change in the voltage magnitude of the infinite bus
"""

mutable struct SourceBusVoltageChange <: Perturbation
    time::Float64
    device::PSY.Source
    signal_index::Int
    ref_value::Float64
end

function get_affect(system::PSY.System, pert::SourceBusVoltageChange)
    device = PSY.get_component(PSY.Source, system, PSY.get_name(pert.device))
    pert.device = device
    return (integrator) -> begin
        control_ref = PSY.get_ext(device)[CONTROL_REFS]
        return control_ref[pert.signal_index] = pert.ref_value
    end
end
