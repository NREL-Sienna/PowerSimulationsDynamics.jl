struct ConstantCurrent end
struct ConstantImpedance end
struct ConstantPower end

const LOAD_MAP = Dict(
    PSY.LoadModels.ConstantCurrent => ConstantCurrent,
    PSY.LoadModels.ConstantImpedance => ConstantImpedance,
    PSY.LoadModels.ConstantPower => ConstantPower,
)
