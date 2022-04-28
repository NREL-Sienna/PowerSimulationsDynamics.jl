import PSSEInterface
reload(PSSEInterface)

powerflow_file='WECC240_v04_psid_comparison.raw'         
dynamic_data_file='WECC240_dynamics_UPV_v04_psid.dyr'         
print('\nInput Files ...')
print(powerflow_file)    
print(dynamic_data_file)

output_dir = '.\\240Bus'        

PSSEInterface.initialize_system(powerflow_file)
# loss_line = 4202-4204

all_perturbations = PSSEInterface.get_line_trips()
#perturbation = ("test", PSSEInterface.LineTrip(2407, 2408, r"""1"""))

PSSEInterface.run_perturbations(all_perturbations, output_dir, powerflow_file, dynamic_data_file, signals = [PSSEInterface.PELEC, PSSEInterface.QELEC])

all_perturbations = PSSEInterface.get_generator_trips()
PSSEInterface.run_perturbations(all_perturbations, output_dir, powerflow_file, dynamic_data_file, signals = [PSSEInterface.PELEC, PSSEInterface.QELEC, PSSEInterface.SPEED])
