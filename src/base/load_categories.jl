abstract type LoadCategory end
struct ConstantCurrent <: LoadCategory end
struct ConstantImpedance <: LoadCategory end
struct ConstantPower <: LoadCategory end

const LOAD_MAP = Dict(
    PSY.LoadModels.ConstantCurrent => ConstantCurrent,
    PSY.LoadModels.ConstantImpedance => ConstantImpedance,
    PSY.LoadModels.ConstantPower => ConstantPower,
)
