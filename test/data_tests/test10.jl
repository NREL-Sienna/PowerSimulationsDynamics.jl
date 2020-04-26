using PowerSystems
using NLsolve
const PSY = PowerSystems

############### Data Network ########################

buses_case10 = buses_multimachine()

branch_case10 = branch_multimachine(buses_case10)

#Trip of Line 1.
branch_case10_fault = branch_multimachine_fault(buses_case10)

loads_case10 = loads_multimachine(buses_case10)


function dyn_gen_second_order(generator)
    return PSY.DynamicGenerator(
        1, #Number
        "Case2_$(get_name(generator))",
        get_bus(generator), #bus
        1.0, # ω_ref,
        1.0, #V_ref
        get_activepower(generator), #P_ref
        get_reactivepower(generator), #Q_ref
        machine_4th(), #machine
        shaft_damping(), #shaft
        avr_type1(), #avr
        tg_none(), #tg
        pss_none(),
    ) #pss
end

function dyn_gen_case10_2(buses)
    return PSY.DynamicGenerator(
        1, #Number
        "Case10Gen2",
        buses[3], #bus
        1.0, # ω_ref,
        1.0, #V_ref
        0.4, #P_ref
        0.0, #Q_ref
        machine_multi(), #machine
        shaft_damping(), #shaft
        avr_none(), #avr
        tg_type2(), #tg
        pss_none(),
    ) #pss
end

#Compute Y_bus after fault
Ybus_fault = Ybus(branch_case10_fault, buses_case10)[:, :]
