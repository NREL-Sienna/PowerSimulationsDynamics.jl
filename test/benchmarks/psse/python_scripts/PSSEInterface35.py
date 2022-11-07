#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
import os
import sys
import re
import matplotlib.pyplot as plt
#
psse_path = r'C:\Program Files\PTI\PSSE35\35.3\PSSBIN'
os.environ['PATH']=os.environ['PATH']+';'+psse_path
#
pssepythonpath = r'C:\Program Files\PTI\PSSE35\35.3\PSSPY39'
sys.path.append(pssepythonpath)
#
import dyntools
import psse35
from psspy import *
from redirect import *
import numpy as np

_i = getdefaultint()
_f = getdefaultreal()

ANGLE = 1
PELEC = 2
QELEC = 3
ETERM = 4
EFD = 5
PMECH = 6
SPEED = 7
XADIFID = 8

csv_file='results.csv'

def initialize_system(powerflow_file, bus_count = 50000):
    psse2py()
    ierr = psseinit(bus_count)
    if ierr != 0:
        raise(Exception("Initialization failed {}".format(ierr)))
    set_NaN_python()
    ierr = read(0, powerflow_file)
    if ierr != 0:
        raise(Exception("Reading File {} failed".format(powerflow_file)))
    ierr = base_frequency(60)
    if ierr != 0:
        raise(Exception("Setting Frequency failed"))

def get_slack_bus():
    ierr, iarray_number = abusint(-1, 2, "NUMBER")
    ierr, iarray = abusint(-1, 2, "TYPE")
    slack_bus = iarray_number[0][iarray[0].index(3)]
    return slack_bus

def convert_power_flow():
    err = 0
    not_finished = True
    while err == 0 and not_finished:
        print('Init Power Flow\n')
        for i in [1,2,3]:
            err = fnsl([0,0,0,1,1,0,0,0])
            err = fnsl([0,0,0,1,1,0,0,0])
            err = fnsl([0,0,0,1,1,0,0,0])
        if err != 0:
            raise(Exception("Power Flow solve failed"))

        print('Load and Gen Convertion to constant impedance\n')
        err = cong(0)
        for i in [1, 2, 3]:
            err = conl(0,1,i,[0,0],[0.0,100.0,0.0,100.0])

        print('Matrix Factorize and Current Power Flow solve\n')
        err = fact()
        err = tysl(0)
        not_finished = False
    if err != 0:
        raise(Exception("Conversion failed"))

def setup_dynamic_simulation(dynamic_data_file, convergence_tolerance = 0.000100, delta_t = 0.008333, network_solve_max_iters = 25):
    print('Load Dynamic data\n')
    err = dyre_new([1,1,1,1], dynamic_data_file,"","","")
    #
    print('Setup Simulation\n')

    err = dynamics_solution_param_2([network_solve_max_iters, _i,_i,_i,_i,_i,_i,_i],[_f, convergence_tolerance, delta_t, _f,_f,_f,_f])
    if err != 0:
        raise(Exception("Dynamics failed to set-up models"))

def setup_output_channels(output, run_name, signals, slack_bus, channel_file='channels.out'):
    print('Set-up output Channels\n')
    results_path = os.path.join(output, run_name)
    channel_file = os.path.join(results_path, channel_file)

    ierr, machine_buses_ix = amachint(-1, 1, "NUMBER")
    ierr, machine_ids = amachchar(-1, 1, "ID")

    ierr = delete_all_plot_channels()
    if ierr != 0:
        raise(Exception("Failed to remove all output channels"))

    channel_ix = 1
    slack_channel = -1
    for (ix, i) in enumerate(machine_buses_ix[0]):
        if i == slack_bus:
            slack_channel = 1*channel_ix
            ierr=machine_array_channel([channel_ix,1,i],machine_ids[0][ix],"")
            if ierr != 0:
                raise(Exception("Channel not set for bus {}, return {}".format(i, ierr)))
        for j in signals:
            if i == slack_bus & j == 1:
                channel_ix += 1
                continue
            ierr=machine_array_channel([channel_ix,j,i],machine_ids[0][ix],"")
            if ierr != 0:
                raise(Exception("Channel not set for bus {}, return {}".format(i, ierr)))
            channel_ix += 1

    if slack_channel < 0:
        raise(Exception("Error determining slack channel"))

    ierr = chsb(0, 1, [-1, -1, -1, 1, 14])
    if ierr != 0:
        raise(Exception("Channels not set for bus quantities"))

    ierr = chsb(0, 1, [-1, -1, -1, 1, 21])
    if ierr != 0:
        raise(Exception("Channels not set for bus quantities"))

    ierr, branch_array = abrnint(-1, _i, _i, 1, 1, ['FROMNUMBER', 'TONUMBER'])
    ierr, branch_ids = abrnchar(-1, _i, _i, 1, 1, "ID")

    for (ix, id) in enumerate(branch_ids[0]):
        bus_from = branch_array[0][ix]
        bus_to = branch_array[1][ix]
        ierr = branch_p_and_q_channel([-1, -1, -1, bus_from, bus_to], id, )
        if ierr != 0:
            raise(Exception("Line channel can't be created {}".format(ierr)))
        ierr = branch_p_and_q_channel([-1, -1, -1, bus_to, bus_from], id, )
        if ierr != 0:
            raise(Exception("Line channel can't be created {}".format(ierr)))

    print('Initialize Simulation\n')
    ierr = strt_2([0, 0], channel_file)
    if ierr != 0:
        raise(Exception("Simulation Initialization failed {}".format(ierr)))

    ierr = set_chnfil_type(0)
    if ierr != 0:
        raise(Exception("Channel file failed to set {}".format(ierr)))

    return slack_channel

def LineTrip(from_node, to_node, id):
    def line_trip():
        ierr = dist_branch_trip(from_node, to_node, id)
        if ierr != 0:
            raise(Exception("Fault can't be created {}".format(ierr)))

    return line_trip

def get_line_trips():
    line_perturbations = []
    ierr, branch_array = abrnint(-1, _i, _i, 1, 1, ['FROMNUMBER', 'TONUMBER'])
    ierr, branch_ids = abrnchar(-1, _i, _i, 1, 1, "ID")
    if ierr != 0:
        raise(Exception("Line data can't be read"))
    for i in range(0, len(branch_array[0])):
        bus_from = branch_array[0][i]
        bus_to = branch_array[1][i]
        line_perturbations.append(("line_{f}_{t}-{id}".format(f = bus_from, t = bus_to, id = branch_ids[0][i]).strip(), LineTrip(bus_from, bus_to, branch_ids[0][i])))
    return line_perturbations

def GenTrip(node, id):
    def gen_trip():
        ierr = dist_machine_trip(node, id)
        if ierr != 0:
            raise(Exception("Fault can't be created {}".format(ierr)))

    return gen_trip

def get_generator_trips():
    gen_perturbations = []
    ierr, machine_buses_ix = amachint(-1, 1, "NUMBER")
    ierr, machine_ids = amachchar(-1, 1, "ID")
    ierr, machine_names = amachchar(-1, 1, "NAME")
    if ierr != 0:
        raise(Exception("Machine data can't be read"))
    for (ix, i) in enumerate(machine_buses_ix[0]):
        gen_perturbations.append(("gen_{}".format(machine_names[0][ix].strip()), GenTrip(i, machine_ids[0][ix])))
    return gen_perturbations


def run_simulation(fault, tspan = (0.0, 10.0), fault_time = 1.0):
    if okstrt()!=0:
        raise(Exception("Simulation Initialization failed"))

    if fault_time > tspan[1]:
        raise(Exception("Fault Time Larger than end of tspan"))

    print('Run Simulation\n')

    ierr = run(0, tspan[0], 1000, 1, 1)
    if ierr != 0:
        raise(Exception("simulation failed"))

    ierr = run(0, fault_time, 999, 1, 1)
    if ierr != 0:
        raise(Exception("simulation failed"))

    fault()

    ierr = run(0, tspan[1], 10000, 1, 1)
    if ierr != 0:
        raise(Exception("simulation failed"))



def process_results(output, run_name, csv_file, slack_channel, make_plots = True, channel_file='channels.out'):
    set_NaN_python()
    results_path = os.path.join(output, run_name)
    channel_file_path = os.path.join(results_path, channel_file)
    print('\nProcessing outputs ...')
    print("\nReading Channel File {}".format(channel_file_path))
    chnfobj = dyntools.CHNF(channel_file_path, outvrsn = 0)
    print('\nChannel File Opened ...')
    csv_file_name = os.path.join(results_path, csv_file)
    slack_channel_id_file_name = os.path.join(results_path, "slack_id.txt")
    chnfobj.csvout(csvfile = csv_file_name)
    with open(slack_channel_id_file_name, 'w') as f:
        f.write("{}".format(slack_channel))
    if make_plots:
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
        if make_plots:
            print("Plotting for {}".format(run_name))
        for k, title in chanid.items():
            if k == 'time':
                continue
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

def close_session():
    ierr = pssehalt_2()
    if ierr != 0:
        raise(Exception("PSSe session not closed succesfully"))

def setup_file_system(output, run_name):
    results_path = os.path.join(output, run_name)
    if (os.path.exists(results_path)):
        print("Overwriting {}".format(results_path))
    else:
        os.makedirs(results_path)

default_signals = [
                ANGLE,
                PELEC,
                QELEC,
                ETERM,
                EFD,
                PMECH,
                SPEED,
                XADIFID
                    ]

def run_perturbations(perturbations, output, powerflow_file, dynamic_data_file, signals = default_signals, convergence_tolerance = 0.001, delta_t = 0.005, make_plots=False):
    progress_output(6,'',[])
    prompt_output(6, '', [])
    if type(perturbations) is tuple:
        perturbations = [perturbations]
    for p in perturbations:
        pertubation_name = p[0]
        print("Simulating {}".format(pertubation_name))
        perturbation_function = p[1]
        setup_file_system(output, pertubation_name)
        initialize_system(powerflow_file)
        slack_bus = get_slack_bus()
        convert_power_flow()
        setup_dynamic_simulation(dynamic_data_file, convergence_tolerance = convergence_tolerance, delta_t = delta_t)
        slack_channel = setup_output_channels(output, pertubation_name, signals, slack_bus)
        run_simulation(perturbation_function, tspan = (0.0, 20.0), fault_time = 1.0)
        process_results(output, pertubation_name, csv_file, slack_channel, make_plots)
    close_session()
