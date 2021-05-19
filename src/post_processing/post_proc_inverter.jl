"""
Function to obtain the output current time series of a Dynamic Inverter model out of the DAE Solution. It receives the simulation inputs,
the dynamic device and bus voltage. It is dispatched for device type to compute the specific current.

"""
function compute_output_current(
    sim::Simulation,
    dynamic_device::G,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
) where {G <: PSY.DynamicInverter}

    #Obtain Data
    sys = get_system(sim.simulation_inputs)

    #Get machine
    filt = PSY.get_filter(dynamic_device)
    Sbase = PSY.get_base_power(sys)
    basepower = PSY.get_base_power(dynamic_device)
    base_power_ratio = basepower / Sbase
    return _filter_current(
        filt,
        PSY.get_name(dynamic_device),
        V_R,
        V_I,
        base_power_ratio,
        sim,
    )
end

"""
Function to obtain the output current time series of a LCL Filter model out of the DAE Solution. It is dispatched via the Filter type.

"""
function _filter_current(
    filt::PSY.LCLFilter,
    name::String,
    V_R::Vector{Float64},
    V_I::Vector{Float64},
    base_power_ratio::Float64,
    sim::Simulation,
)
    ir_filter = post_proc_state_series(sim, (name, :ir_filter))
    ii_filter = post_proc_state_series(sim, (name, :ii_filter))

    return base_power_ratio * ir_filter, base_power_ratio * ii_filter
end
