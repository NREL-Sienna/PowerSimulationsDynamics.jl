OMIB_nodes = nodes_OMIB()

@testset "Dynamic Machines" begin
    Basic = machine_OMIB()
    @test Basic isa PSY.DynamicComponent

    oneDoneQ = machine_4th()
    @test oneDoneQ isa PSY.DynamicComponent

    AndersonFouad = machine_anderson()
    @test AndersonFouad isa PSY.DynamicComponent

    KundurMachine = machine_kundur()
    @test KundurMachine isa PSY.DynamicComponent

    KundurFullMachine = machine_full_kundur()
    @test KundurFullMachine isa PSY.DynamicComponent

    Mach2_benchmark = machine_4th()
    @test Mach2_benchmark isa PSY.DynamicComponent
end

################ Shaft Data #####################
@testset "Dynamic Shaft" begin
    BaseShaft = shaft_damping()
    @test BaseShaft isa PSY.DynamicComponent

    FiveShaft = shaft_fivemass()
    @test FiveShaft isa PSY.DynamicComponent
end

################# PSS Data #####################
@testset "Dynamic PSS" begin
    no_pss = pss_none()
    @test no_pss isa PSY.DynamicComponent
end
################ TG Data #####################
@testset "Dynamic Turbine Governor Constructors" begin
    fixed_tg = tg_none()
    @test fixed_tg isa PSY.DynamicComponent

    typeI_tg = tg_type1()
    @test typeI_tg isa PSY.DynamicComponent

    typeII_tg = tg_type2()
    @test typeII_tg isa PSY.DynamicComponent
end
################ AVR Data #####################
@testset "Dynamic AVR Constructors" begin
    proportional_avr = avr_propr()
    @test proportional_avr isa PSY.DynamicComponent

    fixed_avr = avr_fixed()
    @test fixed_avr isa PSY.DynamicComponent

    typeI_avr = avr_type1()
    @test typeI_avr isa PSY.DynamicComponent

    gen2_avr_benchmark = avr_type2()
    @test gen2_avr_benchmark isa PSY.DynamicComponent
end
######################### Generators ########################
@testset "Dynamic Generators" begin
    #Components for the test
    Basic = machine_OMIB()

    BaseShaft = shaft_damping()

    fixed_avr = avr_fixed()

    proportional_avr = avr_propr()

    fixed_tg = tg_none()

    no_pss = pss_none()

    oneDoneQ = machine_4th()

    Gen1AVR = PSY.DynamicGenerator(
        1, #Number
        "TestGen",
        OMIB_nodes[2],#bus
        1.0, # ω_ref,
        1.05,
        0.4,
        0.0,
        Basic,
        BaseShaft,
        proportional_avr, #avr
        fixed_tg, #tg
        no_pss,
    )
    @test Gen1AVR isa PowerSystems.Component

    Gen1AVRnoAVR = PSY.DynamicGenerator(
        1, #Number
        "TestGen",
        OMIB_nodes[2],#bus
        1.0, # ω_ref,
        1.05,
        0.4,
        0.0,
        Basic,
        BaseShaft,
        fixed_avr, #avr
        fixed_tg, #tg
        no_pss,
    )
    @test Gen1AVRnoAVR isa PowerSystems.Component

    Gen2AVRnoAVR = PSY.DynamicGenerator(
        1, #Number
        "TestGen",
        OMIB_nodes[2],#bus
        1.0, # ω_ref,
        1.02,
        0.4,
        0.0,
        oneDoneQ,
        BaseShaft,
        fixed_avr, #avr
        fixed_tg, #tg
        no_pss,
    )
    @test Gen2AVRnoAVR isa PowerSystems.Component

    Gen2AVR = PSY.DynamicGenerator(
        1, #Number
        "TestGen",
        OMIB_nodes[2],#bus
        1.0, # ω_ref,
        1.02,
        0.4,
        0.0,
        oneDoneQ,
        BaseShaft,
        proportional_avr, #avr
        fixed_tg, #tg
        no_pss,
    )
    @test Gen2AVR isa PowerSystems.Component

end
