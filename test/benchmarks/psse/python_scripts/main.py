import PSSEInterface
import os

powerflow_file='ThreeBusRenewable.raw'         
dynamic_data_file='ThreeBus_REN_A_NOFREQFLAG_with_REF_FLAG.dyr'      
print '\nInput Files ...'
print powerflow_file    
print dynamic_data_file

channel_file='channels.out'  
csv_file='results.csv'
output = '.'        
print '\Output files...'
print channel_file
print csv_file

perturbations = [("line101_103", PSSEInterface.LineFault(101, 103))]

PSSEInterface.init_system(powerflow_file)
slack_bus = PSSEInterface.get_slack_bus()

for p in perturbations:
    results_path = os.path.join(output, p[0])
    os.mkdir(results_path)
    channel_file_path = os.path.join(results_path, channel_file)
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
    slack_channel = PSSEInterface.SetUpOutputChannels(channel_file_path, signals, slack_bus)
    PSSEInterface.RunSimulation(p[1], tspan = (0.0, 10.0), fault_time = 1.0)
    PSSEInterface.ProcessResults(results_path, channel_file, csv_file, slack_channel)
    PSSEInterface.CloseSession()