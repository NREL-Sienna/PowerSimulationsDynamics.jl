isdefined(Base, :__precompile__) && __precompile__()
module PowerSimulationsDynamics

####################################### Structs Exports ####################################

# Base Exports
export Simulation
export Simulation!
export execute!
export NetworkSwitch
export ControlReferenceChange
export BranchTrip
# export BusTrip

export ImplicitModel
export MassMatrixModel

# Export for routines
export small_signal_analysis
export get_state_series
export get_voltagemag_series
export print_device_states
export get_initial_conditions
export get_activepower_series
export get_reactivepower_series

####################################### Package Imports ####################################
import Logging
import InfrastructureSystems
import SciMLBase
import DiffEqBase
import ForwardDiff
import SparseArrays: SparseMatrixCSC, sparse
import LinearAlgebra
import Base.to_index
import NLsolve
import Base.ImmutableDict
import PowerSystems
const PSY = PowerSystems
const IS = InfrastructureSystems
const PSID = PowerSimulationsDynamics

using DocStringExtensions

#@template (FUNCTIONS, METHODS) = """
#                                 $(TYPEDSIGNATURES)
#                                 $(DOCSTRING)
#                                 """

#Structs for General Devices and System
include("base/definitions.jl")
include("base/ports.jl")
include("base/perturbations.jl")
include("base/file_system.jl")
include("base/simulation_model.jl")
include("base/simulation_inputs.jl")
include("base/device_indexing.jl")
include("base/mass_matrix.jl")
include("base/frequency_reference.jl")
include("base/simulation.jl")
include("base/supplemental_accesors.jl")
include("base/simulation_initialization.jl")
include("base/small_signal.jl")

#Common Models
include("models/branch.jl")
include("models/device.jl")
include("models/kirchoff_laws.jl")
include("models/dynline_model.jl")
include("models/ref_transformations.jl")

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

#Utils
include("utils/plot_utils.jl")
include("utils/immutable_dicts.jl")
include("utils/print.jl")
include("utils/kwargs_check.jl")
include("utils/logging.jl")

end # module
