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

function _record_change!(
    ybus::SparseArrays.SparseMatrixCSC{Float64, Int},
    bus_from_no::Int,
    bus_to_no::Int,
    n_buses::Int,
    Y11_real::Float64,
    Y11_imag::Float64,
    Y22_real::Float64,
    Y22_imag::Float64,
    Y12_real::Float64,
    Y21_real::Float64,
    Y12_imag::Float64,
    Y21_imag::Float64,
)

    # First Quadrant Real Part Changes
    ybus[bus_from_no, bus_from_no] += Y11_real
    ybus[bus_from_no, bus_to_no] += Y12_real
    ybus[bus_to_no, bus_from_no] += Y21_real
    ybus[bus_to_no, bus_to_no] += Y22_real

    # Second Quadrant Imag Part changes
    ybus[bus_from_no, bus_from_no + n_buses] += Y11_imag
    ybus[bus_from_no, bus_to_no + n_buses] += Y12_imag
    ybus[bus_to_no, bus_from_no + n_buses] += Y21_imag
    ybus[bus_to_no, bus_to_no + n_buses] += Y22_imag

    # Third Quadrant -1*Imag Part changes
    ybus[bus_from_no + n_buses, bus_from_no] -= Y11_imag
    ybus[bus_from_no + n_buses, bus_to_no] -= Y12_imag
    ybus[bus_to_no + n_buses, bus_from_no] -= Y21_imag
    ybus[bus_to_no + n_buses, bus_to_no] -= Y22_imag

    # Fourth Quadrant Real Part Changes
    ybus[bus_from_no + n_buses, bus_from_no + n_buses] += Y11_real
    ybus[bus_from_no + n_buses, bus_to_no + n_buses] += Y12_real
    ybus[bus_to_no + n_buses, bus_from_no + n_buses] += Y21_real
    ybus[bus_to_no + n_buses, bus_to_no + n_buses] += Y22_real
    return
end

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{Float64, Int},
    b::PSY.ACBranch,
    num_bus::Dict{Int, Int},
    mult::Float64,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[PSY.get_number(arc.from)]
    bus_to_no = num_bus[arc.to.number]
    n_buses = length(num_bus)

    Y_l = (1 / (PSY.get_r(b) + PSY.get_x(b) * 1im))

    Y11_real = mult * real(Y_l)
    Y11_imag = mult * (imag(Y_l) + PSY.get_b(b).from)
    Y22_real = mult * real(Y_l)
    Y22_imag = mult * (imag(Y_l) + PSY.get_b(b).to)
    Y12_real = Y21_real = -mult * real(Y_l)
    Y12_imag = Y21_imag = -mult * imag(Y_l)

    _record_change!(
        ybus,
        bus_from_no,
        bus_to_no,
        n_buses,
        Y11_real,
        Y11_imag,
        Y22_real,
        Y22_imag,
        Y12_real,
        Y21_real,
        Y12_imag,
        Y21_imag,
    )
    return
end

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{Float64, Int},
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
    Y21 = Y12 = (-Y_t * c)
    Y22 = Y_t + (1im * b)

    Y11_real = mult * real(Y11)
    Y11_imag = mult * imag(Y12)
    Y22_real = mult * real(Y22)
    Y22_imag = mult * imag(Y22)
    Y12_real = -mult * real(Y12)
    Y21_real = -mult * real(Y21)
    Y12_imag = -mult * imag(Y12)
    Y21_imag = -mult * imag(Y21)

    _record_change!(
        ybus,
        bus_from_no,
        bus_to_no,
        n_buses,
        Y11_real,
        Y11_imag,
        Y22_real,
        Y22_imag,
        Y12_real,
        Y21_real,
        Y12_imag,
        Y21_imag,
    )
    return
end

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{Float64, Int},
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
    Y12 = (-Y_t / c_tap)
    Y21 = (-Y_t / tap)
    Y22 = Y_t

    Y11_real = mult * real(Y11)
    Y11_imag = mult * imag(Y12)
    Y22_real = mult * real(Y22)
    Y22_imag = mult * imag(Y22)
    Y12_real = -mult * real(Y12)
    Y21_real = -mult * real(Y21)
    Y12_imag = -mult * imag(Y12)
    Y21_imag = -mult * imag(Y21)

    _record_change!(
        ybus,
        bus_from_no,
        bus_to_no,
        n_buses,
        Y11_real,
        Y11_imag,
        Y22_real,
        Y22_imag,
        Y12_real,
        Y21_real,
        Y12_imag,
        Y21_imag,
    )

    return
end

function ybus_update!(
    ybus::SparseArrays.SparseMatrixCSC{Float64, Int},
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
    ybus_rectangular::SparseArrays.SparseMatrixCSC{Float64, Int}
    function NetworkSwitch(
        time::Float64,
        ybus::SparseArrays.SparseMatrixCSC{Complex{Float64}, Int},
    )
        n_bus = size(ybus)[1]
        if n_bus < 15_000
            I = PSY._goderya(ybus)
            if length(Set(I)) != n_bus
                error("The Ybus provided has islands and can't be used in the simulation")
            end
        else
            @warn(
                "The number of buses in the Ybus provided is $n_bus and is too large to be verified for connectivity"
            )
        end
        # TODO: Improve performance here
        ybus_rect = transform_ybus_to_rectangular(ybus)
        new(time, ybus_rect)
    end
end

function NetworkSwitch(time::Float64, ybus::PSY.Ybus)
    return NetworkSwitch(time, ybus.data)
end

function get_affect(::PSY.System, pert::NetworkSwitch)
    return (integrator) -> begin
        # TODO: This code can be more performant using SparseMatrix methods
        for (i, v) in enumerate(pert.ybus_rectangular)
            integrator.p.ybus_reactangular[i] = v
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
