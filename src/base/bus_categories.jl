abstract type BusCategory end
struct PVBus <: BusCategory end
struct PQBus <: BusCategory end
struct SLACKBus <: BusCategory end

const BUS_MAP = Dict(
    PSY.BusTypes.REF => SLACKBus,
    PSY.BusTypes.PV => PVBus,
    PSY.BusTypes.SLACK => SLACKBus,
    PSY.BusTypes.PQ => PQBus,
)
