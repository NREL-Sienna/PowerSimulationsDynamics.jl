import PSSEInterface
reload(PSSEInterface)

powerflow_file='WECC240_v04_DPV_RE20_v33_6302_xfmr_DPbuscode_PFadjusted_V32.raw'         
dynamic_data_file='WECC240_dynamics_UPV_v04_all.dyr'         
print '\nInput Files ...'
print powerflow_file    
print dynamic_data_file

output_dir = '.\\240Bus'        

PSSEInterface.initialize_system(powerflow_file)
all_perturbations = PSSEInterface.get_generator_trips()[0:2]
#perturbation = ("Line", PSSEInterface.LineTrip(2407, 2408, r"""1"""))

PSSEInterface.run_perturbations(all_perturbations, output_dir, powerflow_file, dynamic_data_file, signals = [PSSEInterface.PELEC, PSSEInterface.QELEC])