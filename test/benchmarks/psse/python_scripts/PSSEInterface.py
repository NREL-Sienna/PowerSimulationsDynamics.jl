#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
import os
import sys
import re
import matplotlib.pyplot as plt
#
psse_path = r'C:\Program Files (x86)\PTI\PSSE34\PSSBIN'
os.environ['PATH']=os.environ['PATH']+';'+psse_path
#
pssepythonpath = r'C:\Program Files (x86)\PTI\PSSE34\PSSPY27'
sys.path.append(pssepythonpath)
#
import dyntools
from psspy import *
from redirect import *
import numpy as np

_i = getdefaultint()  
_f = getdefaultreal()   

def init_system(powerflow_file):  
    psse2py()              
    psseinit(50)            
    read(0, powerflow_file)                         
    ierr = base_frequency(60) 

def get_slack_bus():
    ierr, iarray_number = abusint(-1, 2, "NUMBER") 
    ierr, iarray = abusint(-1, 2, "TYPE")
    slack_bus = iarray_number[0][iarray[0].index(3)]
    return slack_bus

def PowerFlowConvert():
    err = 0
    not_finished = True
    while err == 0 and not_finished:
        print '\nInit Power Flow ...' 
        for i in [1,2,3]:  
            err = fnsl([0,0,0,1,1,0,0,0])     
        
        print '\nLoad and Gen Convertion ...'   
        err = cong(0)
        for i in [1, 2, 3]:
            err = conl(0,1,i,[0,0],[100.0,0.0,0.0,100.0])

        print '\nMatrix Factorize and Current Power Flow solve ...'   
        err = fact()     
        err = tysl(0)
        not_finished = False
    if err != 0:
        raise(Exception("Conversion failed"))

def SetUpDynamicSimulation(dynamic_data_file, convergence_tolerance = 0.0001, delta_t = 0.005, network_solve_max_iters = _i):
    print '\nLoad DynamicDdata ...'   
    err = dyre_new([1,1,1,1], dynamic_data_file,"","","")
    #
    print '\nSetup Simulation ...'   
    
    err = dynamics_solution_param_2([network_solve_max_iters, _i,_i,_i,_i,_i,_i,_i],[_f,_f,delta_t,convergence_tolerance,_f,_f,_f,_f])
    if err != 0:
        raise(Exception("Dynamics failed to set-up models"))

def SetUpOutputChannels(output, run_name, signals, slack_bus, channel_file='channels.out'):
    print '\nSet-up output Channels ...'   

    results_path = os.path.join(output, run_name)
    channel_file = os.path.join(results_path, channel_file)
    
    ierr, machine_buses = agenbusint(-1, 1, "NUMBER")
    ierr = delete_all_plot_channels()
    if ierr != 0:
        raise(Exception("Failed to remove all output channels"))

    channel_ix = 1
    slack_channel = -1
    for i in machine_buses[0]:
        for j in signals:
            ierr=machine_array_channel([channel_ix,j,i],r"""1""","")
            if ierr != 0:
                raise(Exception("Channel not set, return {}".format(ierr)))
            if i == slack_bus and j == 1:
                slack_channel = 1*channel_ix
            channel_ix += 1


    ierr = chsb(0, 1, [-1, -1, -1, 1, 13])
    if ierr != 0:
        raise(Exception("Channels not set for bus quantities"))
    
    print '\nInitialize Simulation ...'   
    ierr = strt_2([0, 0], channel_file)
    if ierr != 0:
        raise(Exception("Simulation Initialization failed"))

    if slack_channel < 0:
        raise(Exception("Error determining slack channel"))

    return slack_channel
       
def LineFault(from_node, to_node):
    def line_trip(): 
        ierr = dist_branch_trip(from_node, to_node, r"""1""")   
        if ierr != 0:
            raise(Exception("Fault can't be created"))   
    
    return line_trip

def get_line_trips():
    line_perturbations = []
    ierr, branch_array = abrnint(-1, _i, _i, 1, 1, ['FROMNUMBER', 'TONUMBER'])
    if ierr != 0:
        raise(Exception("Line data can't be read")) 
    for i in range(0, len(branch_array[0])):
        bus_from = branch_array[0][i]
        bus_to = branch_array[1][i]
        line_perturbations.append(("line_{f}_{t}".format(f = bus_from, t = bus_to), LineFault( bus_from, bus_to)))
    return line_perturbations

def RunSimulation(fault, tspan = (0.0, 10.0), fault_time = 1.0):
    if okstrt()!=0:
        raise(Exception("Simulation Initialization failed"))
    
    print '\nRun Simulation ...'   

    run(0, fault_time, 100000, 1, 1) 
    
    fault()
    
    ierr = run(0, tspan[1], 10000, 1, 1)
    if ierr != 0:
        raise(Exception("simulation failed"))   

def ProcessResults(output, run_name, csv_file, slack_channel, channel_file='channels.out'):
    results_path = os.path.join(output, run_name)
    channel_file_path = os.path.join(results_path, channel_file)
    print '\nProcessing outputs ...'   
    print "\nReading Channel File {}".format(channel_file_path)
    chnfobj = dyntools.CHNF(channel_file_path)
    print '\nChannel File Opened ...'
    _, chanid, chandata = chnfobj.get_data()
    t = chandata['time']
    export = np.array(t)
    header = "t,"
    tspan = (min(t), max(t))
    
    plots_path = os.path.join(results_path, "plots")
    if (os.path.exists(plots_path)):
        print("Overwriting {}".format(plots_path))
    else:    
        os.mkdir(plots_path)

    for k, title in chanid.items():
        if k == 'time':
            continue
        print(k)
        if bool(re.match(r'ANGL', title)):
            result = np.array(chandata[k]) - np.array(chandata[slack_channel])
        else:
            result = np.array(chandata[k])
        export = np.vstack([export, result])
        header = header+"{},".format(title)
        fig = plt.figure()
        fig.patch.set_facecolor('0.8')
        plt.plot(t,result,linestyle='-', linewidth=1, color='red')
        plt.grid(linestyle='--', color='grey',linewidth=0.5)
        plt.xlabel("Time, s")
        plt.ylabel(title)
        plt.legend()
        axes = plt.gca()
        axes.set_facecolor('w')
        axes.set_ylim([min(result),max(result)])
        axes.set_xlim([tspan[0],tspan[1]])
        plt.savefig(os.path.join(plots_path, "{}.pdf".format(title)),dpi=fig.dpi,facecolor='0.8')
        plt.close()
    np.savetxt(os.path.join(results_path, csv_file), np.transpose(export), delimiter=",", fmt="%.8f", header=header, comments="")

def CloseSession():
    ierr = pssehalt_2()
    if ierr != 0:
        raise(Exception("PSSe session not closed succesfully"))

def FileSystemSetUp(output, run_name):
    results_path = os.path.join(output, run_name)
    if (os.path.exists(results_path)):
        print("Overwriting {}".format(results_path))
    else:    
        os.makedirs(results_path)  