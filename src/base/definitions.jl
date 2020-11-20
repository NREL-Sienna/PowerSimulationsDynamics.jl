#=
Inner variables are used for keeping track of internal variables
of different dynamic devices. These are not states, since can be computed
algebraically using states and parameters, and are only used internally in each
component. In such case, it is possible to avoid creating dummy states that
must be added to the vector of states, by handling those as inner variables.
In some cases, depending on the complexity of the model, some inner vars are
also defined as states, but for the sake of generality are also added as inner
variables, since some models may not treat such variables as states.
=#

"""
Generator Inner Vars:
*  τe_var :: Electric torque
*  τm_var :: Mechanical torque
*  Vf_var :: Field voltage
*  V_pss_var :: Additional PSS voltage
*  VR_gen_var :: Real part of the terminal voltage
*  VI_gen_var :: Imaginary part of the terminal voltage
*  ψd_var :: Stator Flux (if defined) in the d-axis
*  ψq_var :: Stator Flux (if defined) in the q-axis
"""
@enum generator_inner_vars begin
    τe_var = 1
    τm_var = 2
    Vf_var = 3
    V_pss_var = 4
    VR_gen_var = 5
    VI_gen_var = 6
    ψd_var = 7
    ψq_var = 8
    Xad_Ifd_var = 9
end

Base.to_index(ix::generator_inner_vars) = Int(ix)

"""
Inverter Inner Vars:
@enum inverter_inner_vars begin
md_var :: Modulation signal on the d-component
mq_var :: Modulation signal on the q-component
Vdc_var :: DC voltage supplied by the DC source
Vr_filter_var :: Voltage seen in the capacitor of the filter in the R-component
Vi_filter_var :: Voltage seen in the capacitor of the filter in the I-component
ω_freq_estimator_var :: Frequency estimated by the frequency estimator.
V_oc_var :: Control voltage supplied from the outer loop control to the inner loop
ω_oc_var :: Control frequency supplied from the outer loop control the inner loop
θ_oc_var :: Variation of the angle (PLL or VSM) of the inverter
VR_inv_var :: Real terminal voltage on the inverter
VI_inv_var :: Imaginary terminal voltage on the inverter
Vr_cnv_var :: Voltage supplied from the converter in the R-component
Vi_cnv_var :: Voltage supplied from the converter in the I-component
"""
@enum inverter_inner_vars begin
    md_var = 1
    mq_var = 2
    Vdc_var = 3
    Vr_filter_var = 4
    Vi_filter_var = 5
    ω_freq_estimator_var = 6
    V_oc_var = 7
    ω_oc_var = 8
    θ_oc_var = 9
    VR_inv_var = 10
    VI_inv_var = 11
    Vr_cnv_var = 12
    Vi_cnv_var = 13
    P_ES_var = 14 #Power from energy source
end

Base.to_index(ix::inverter_inner_vars) = Int(ix)

@enum dq_ref begin
    d = 1
    q = 2
end
@enum RI_ref begin
    R = 1
    I = 2
end

Base.to_index(ix::dq_ref) = Int(ix)
Base.to_index(ix::RI_ref) = Int(ix)

const V_ref_index = 1
const ω_ref_index = 2
const P_ref_index = 3
const Q_ref_index = 4

const MAPPING_DICT = Dict{String, Dict{Symbol, Int64}}

const LOCAL_STATE_MAPPING = "local_state_mapping"
const INPUT_PORT_MAPPING = "input_port_mapping"
const PORTS = "ports"
const INNER_VARS = "inner_vars"
const CONTROL_REFS = "control_refs"

const SIMULATION_ACCEPTED_KWARGS =
    [:initialize_simulation, :system_to_file, :file_level, :console_level]
const SMALL_SIGNAL_ACCEPTED_KWARGS = [:reset_simulation!]
const RELAXED_NL_SOLVE_TOLERANCE = :1e-6
const STRICT_NL_SOLVE_TOLERANCE = :1e-9
const MINIMAL_ACCEPTABLE_NL_SOLVE_TOLERANCE = :1e-3

const SIMULATION_LOG_FILENAME = "power-simulations-dynamics.log"

"""
Defines the status of the simulation object
"""
@enum BUILD_STATUS begin
    BUILT = 0
    BUILD_INCOMPLETE = 1
    BUILD_FAILED = 2
    SIMULATION_STARTED = 3
    SIMULATION_FINALIZED = 4
    SIMULATION_FAILED = 5
    CONVERTED_FOR_SMALL_SIGNAL = 6
end
