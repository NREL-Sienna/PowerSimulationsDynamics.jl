using PowerSystems
const PSY = PowerSystems

############### Data Network ########################
include(joinpath(dirname(@__FILE__), "dynamic_test_data.jl"))
include(joinpath(dirname(@__FILE__), "data_utils.jl"))
############### Data Network ########################
threebus_file_dir = joinpath(dirname(@__FILE__), "ThreeBusInverter.raw")
threebus_sys = System(
    PowerModelsData(threebus_file_dir),
    runchecks = false,
    unit_system = "device_base",
)
add_source_to_ref(threebus_sys)
dyn_branch = DynamicBranch(get_component(Branch, threebus_sys, "3"))
add_component!(threebus_sys, dyn_branch)

############### Data devices ########################

function dyn_gen_second_order(generator)
    return PSY.DynamicGenerator(
        generator,
        1.0, # ω_ref,
        machine_oneDoneQ(), #machine
        shaft_no_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(), #pss
    )
end

function inv_case78(static_device)
    return PSY.DynamicInverter(
        static_device,
        1.0, # ω_ref,
        converter_high_power(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source_lv(), #dc source
        pll(), #pll
        filt(), #filter
    )
end

for g in get_components(Generator, threebus_sys)
    if get_number(get_bus(g)) == 102
        case_gen = dyn_gen_second_order(g)
        add_component!(threebus_sys, case_gen)
    elseif get_number(get_bus(g)) == 103
        case_inv = inv_case78(g)
        add_component!(threebus_sys, case_inv)
    end
end

#Create Ybus_Fault
sys3 = System(PowerModelsData(threebus_file_dir), runchecks = false)
add_source_to_ref(sys3)
remove_component!(Line, sys3, "3")
#Create Ybus_Fault
fault_branches2 = get_components(Line, sys3)
for br in fault_branches2
    if get_name(br) == "1"
        br.r = 3 * br.r
        br.x = 3 * br.x
        b_new = (from = br.b.from / 3, to = br.b.to / 3)
        br.b = b_new
    end
end
Ybus_fault = Ybus(sys3).data;
