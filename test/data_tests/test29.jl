using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusRenewable.raw")

### Case 2 Generators ###

function inv_generic_renewable(static_device, F_Flag)
    if F_Flag == 0
        return PSY.DynamicInverter(
            get_name(static_device),
            1.0, #ω_ref
            converter_regca(), #converter
            outer_control_TypeAB(), #outercontrol
            inner_ctrl_typeB(), #inner_control
            dc_source_lv(),
            no_pll(),
            filt_current(),
        )
    else
        return PSY.DynamicInverter(
            get_name(static_device),
            1.0, #ω_ref
            converter_regca(), #converter
            outer_control_TypeAB_FreqFlag(), #outercontrol
            inner_ctrl_typeB(), #inner_control
            dc_source_lv(),
            no_pll(),
            filt_current(),
        )
    end
end
