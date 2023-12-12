abstract type BusCategory end
struct PVBus <: BusCategory end
struct PQBus <: BusCategory end
struct SLACKBus <: BusCategory end

const BUS_MAP = Dict(
    PSY.ACBusTypes.REF => SLACKBus,
    PSY.ACBusTypes.PV => PVBus,
    PSY.ACBusTypes.SLACK => SLACKBus,
    PSY.ACBusTypes.PQ => PQBus,
)
