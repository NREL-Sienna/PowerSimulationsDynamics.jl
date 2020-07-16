function ybus!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.ACBranch,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]

    Y_l = (1 / (PSY.get_r(b) + PSY.get_x(b) * 1im))
    Y11 = Y_l + (1im * PSY.get_b(b).from)

    ybus[bus_from_no, bus_from_no] += mult * Y11
    Y12 = -Y_l
    ybus[bus_from_no, bus_to_no] += mult * Y12
    #Y21 = Y12
    ybus[bus_to_no, bus_from_no] += mult * Y12
    Y22 = Y_l + (1im * PSY.get_b(b).to)
    ybus[bus_to_no, bus_to_no] += mult * Y22
    return
end

function ybus!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.TapTransformer,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[arc.from.number]
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

function ybus!(
    ybus::SparseMatrixCSC{ComplexF64, Int64},
    b::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int64, Int64},
    mult::Float64 = 1.0,
)
    arc = PSY.get_arc(b)
    bus_from_no = num_bus[arc.from.number]
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

function ybus!(
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
