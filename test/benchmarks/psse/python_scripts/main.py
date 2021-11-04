import PSSEInterface
import os

powerflow_file='C:\Users\jlara\Desktop\case\mod_ren14.raw'         
dynamic_data_file='C:\Users\jlara\Desktop\case\dyn_data_14bus_mod.dyr'         
print '\nInput Files ...'
print powerflow_file    
print dynamic_data_file

csv_file='results.csv'
output = '.\\14-Bus'        
print '\Output files...'
print csv_file

PSSEInterface.init_system(powerflow_file)
perturbations = PSSEInterface.get_line_trips()
slack_bus = PSSEInterface.get_slack_bus()

for p in perturbations:
    PSSEInterface.FileSystemSetUp(output, p[0])
    PSSEInterface.init_system(powerflow_file)
    PSSEInterface.PowerFlowConvert()
    PSSEInterface.SetUpDynamicSimulation(dynamic_data_file, convergence_tolerance = 0.0001, delta_t = 0.005)

    signals = [
            1, #ANGLE,
            2, #PELEC,
            3, #QELEC, 
            4, #ETERM,
            5, #EFD,
            6, #PMECH,
            7, #SPEED,
            #8, #XADIFID, 
        ]
    slack_channel = PSSEInterface.SetUpOutputChannels(output, p[0], signals, slack_bus)
    PSSEInterface.RunSimulation(p[1], tspan = (0.0, 10.0), fault_time = 1.0)
    PSSEInterface.ProcessResults(output, p[0], csv_file, slack_channel)
    PSSEInterface.CloseSession()