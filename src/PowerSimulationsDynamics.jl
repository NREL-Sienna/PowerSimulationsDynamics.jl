isdefined(Base, :__precompile__) && __precompile__()
module PowerSimulationsDynamics

####################################### Structs Exports ####################################

# Base Exports
export ResidualModel
export MassMatrixModel

export Simulation
export Simulation!
export SimulationResults

export execute!
export read_results

# Frequency reference options
export ReferenceBus
export ConstantFrequency

# export perturbations
export NetworkSwitch
export ControlReferenceChange
export BranchTrip
export BranchImpedanceChange
export SourceBusVoltageChange
export GeneratorTrip
export LoadTrip
export LoadChange
export PerturbState
# export BusTrip

# Export for routines
export get_jacobian
export small_signal_analysis
export summary_participation_factors
export summary_eigenvalues
export get_state_series
export get_voltage_magnitude_series
export get_voltage_angle_series
export show_states_initial_value
export read_initial_conditions
export get_real_current_series
export get_imaginary_current_series
export get_activepower_series
export get_reactivepower_series
export get_source_real_current_series
export get_source_imaginary_current_series
export get_field_current_series
export get_field_voltage_series
export get_pss_output_series
export get_mechanical_torque_series
export get_real_current_branch_flow
export get_imaginary_current_branch_flow
export get_activepower_branch_flow
export get_reactivepower_branch_flow
export get_frequency_series
export get_setpoints
export get_solution
export is_valid
export transform_load_to_constant_impedance
export transform_load_to_constant_current
export transform_load_to_constant_power
#export get_parameter_values
#export get_forward_function
#export get_gradient_function

####################################### Package Imports ####################################
import Logging
import InfrastructureSystems
import SciMLBase
import DataStructures: OrderedDict
import DataFrames: DataFrame
import Random
import ForwardDiff
import SparseArrays
import LinearAlgebra
import Base.to_index
import NonlinearSolve
import PrettyTables
import Base.ImmutableDict
import PowerSystems
import PowerFlows
import PowerNetworkMatrices
import TimerOutputs
import FastClosures: @closure
import Enzyme
Enzyme.API.runtimeActivity!(true)
Enzyme.API.looseTypeAnalysis!(true)  #Required for using component arrays with Enzyme
Enzyme.API.maxtypeoffset!(1024)
Enzyme.API.maxtypedepth!(20)
#Enzyme.API.runtimeActivity!(true)
#= Generally, the preferred solution to these type of activity unstable codes should be to make 
your variables all activity-stable (e.g. always containing differentiable memory or always
 containing non-differentiable memory). However, with care, Enzyme does support "Runtime Activity" 
 as a way to differentiate these programs without having to modify your code. =#
#Enzyme.API.looseTypeAnalysis!(true)    #Best guess if it cannot determine the type; this was needed for using ComponentArrays/ODEFunction
#Enzyme.API.maxtypedepth!(20)    
import ChainRulesCore
import ComponentArrays

const PSY = PowerSystems
const IS = InfrastructureSystems
const PSID = PowerSimulationsDynamics
const PF = PowerFlows
const PNM = PowerNetworkMatrices
const CRC = ChainRulesCore

using DocStringExtensions

#@template (FUNCTIONS, METHODS) = """
#                                 $(TYPEDSIGNATURES)
#                                 $(DOCSTRING)
#                                 """

#Structs for General Devices and System
include("base/definitions.jl")
include("base/ports.jl")
include("base/auto_abstol.jl")
include("base/bus_categories.jl")
include("base/device_wrapper.jl")
include("base/branch_wrapper.jl")
include("base/frequency_reference.jl")
include("base/simulation_model.jl")
include("base/simulation_inputs.jl")
include("base/perturbations.jl")
include("base/caches.jl")
include("base/system_model.jl")
include("base/jacobian.jl")
include("base/mass_matrix.jl")
include("base/simulation_results.jl")
include("base/simulation.jl")
include("base/global_variables.jl")
include("base/supplemental_accesors.jl")
include("base/nlsolve_wrapper.jl")
include("base/simulation_initialization.jl")
include("base/small_signal.jl")
include("base/model_validation.jl")
#include("base/sensitivity_analysis.jl")

#Common Models
include("models/branch.jl")
include("models/device.jl")
include("models/network_model.jl")
include("models/dynline_model.jl")
include("models/ref_transformations.jl")
include("models/common_controls.jl")

#Generator Component Models
include("models/generator_models/machine_models.jl")
include("models/generator_models/pss_models.jl")
include("models/generator_models/avr_models.jl")
include("models/generator_models/tg_models.jl")
include("models/generator_models/shaft_models.jl")

#Inverter Component Models
include("models/inverter_models/DCside_models.jl")
include("models/inverter_models/filter_models.jl")
include("models/inverter_models/frequency_estimator_models.jl")
include("models/inverter_models/outer_control_models.jl")
include("models/inverter_models/inner_control_models.jl")
include("models/inverter_models/converter_models.jl")

#Injection Models
include("models/load_models.jl")
include("models/source_models.jl")

#Saturation Models
include("models/saturation_models.jl")

#Initialization Parameters
include("initialization/init_device.jl")

#Initialization Generator
include("initialization/generator_components/init_machine.jl")
include("initialization/generator_components/init_shaft.jl")
include("initialization/generator_components/init_avr.jl")
include("initialization/generator_components/init_tg.jl")
include("initialization/generator_components/init_pss.jl")

#Initialization Inverter
include("initialization/inverter_components/init_filter.jl")
include("initialization/inverter_components/init_DCside.jl")
include("initialization/inverter_components/init_converter.jl")
include("initialization/inverter_components/init_frequency_estimator.jl")
include("initialization/inverter_components/init_inner.jl")
include("initialization/inverter_components/init_outer.jl")

#System Model
include("models/system.jl")

#PostProcessing
include("post_processing/post_proc_common.jl")
include("post_processing/post_proc_generator.jl")
include("post_processing/post_proc_inverter.jl")
include("post_processing/post_proc_results.jl")
include("post_processing/post_proc_loads.jl")
include("post_processing/post_proc_source.jl")

#Utils
include("utils/psy_utils.jl")
include("utils/immutable_dicts.jl")
include("utils/print.jl")
include("utils/kwargs_check.jl")
include("utils/logging.jl")
include("utils/parameters.jl")

end # module
