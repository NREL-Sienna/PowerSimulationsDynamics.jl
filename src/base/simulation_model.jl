abstract type ModelVariations end

abstract type SimulationModel <: ModelVariations end
struct MassMatrixModel <: SimulationModel end
struct ResidualModel <: SimulationModel end

abstract type DelayModel <: ModelVariations end
struct HasDelays <: DelayModel end
struct NoDelays <: DelayModel end
