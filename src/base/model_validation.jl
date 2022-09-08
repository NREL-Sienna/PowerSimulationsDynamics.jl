# By default all combinations are valid, if a combination of parameters can't be modeled
# needs to be added here to dispatch on the invalid dynamic component types.

is_valid(::PSY.DynamicInverter) = nothing
is_valid(::PSY.DynamicGenerator) = nothing

function is_valid(
    device::PSY.DynamicGenerator{PSY.BaseMachine, S, PSY.AVRFixed, TG, P},
) where {S <: PSY.Shaft, TG <: PSY.TurbineGov, P <: PSY.PSS}
    error("Device $(PSY.get_name(device)) uses a BaseMachine model with an AVR")
    return
end

function is_valid(
    device::PSY.DynamicInverter{
        PSY.AverageConverter,
        PSY.OuterControl{PSY.ActivePowerPI, PSY.ReactivePowerPI},
        PSY.VoltageModeControl,
        DC,
        P,
        F,
    },
) where {DC <: PSY.DCSource, P <: PSY.FrequencyEstimator, F <: PSY.Filter}
    error(
        "Device $(PSY.get_name(device)) uses a ActivePowerPI and ReactivePowerPI model with a VoltageModeControl",
    )
    return
end
