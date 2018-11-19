#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: wtp
"""

"""
structure:
   1) in first timestep: read control file
   2) re-write input files for storage simulation
   3) run storage simulation
   4) read results and gather well data
   5) return one pressure and the storage flowrate to tespy

"""

class ctrl_sto:
    
    def __init__(self):
        self.path = ''
        self.simulator = ''
        self.simulator_path = ''
        self.simulation_title = ''
        self.restart_id = ''
        self.well_names = []
        self.well_lower_BHP = []
        self.well_upper_BHP = []
        self.keep_ecl_logs = False
    
    def set_ecl_log_flag(self, a_flag):
        self.keep_ecl_logs = a_flag

    def set_path(self, a_path):
        self.path = a_path
    
    def set_simulator(self, simulator):
        self.simulator = simulator
    
    def set_restart_id(self, a_number):
        self.restart_id = str(a_number)

    def set_simulation_title(self, title):
        self.simulation_title = title

    def set_simulator_path(self, simulator_path):
        self.simulator_path = simulator_path

    def add_well_name(self, a_name):
        self.well_names.append(a_name)
    
    def add_well_lower_BHP_val(self, value):
        self.well_lower_BHP.append(value)

    def add_well_upper_BHP_val(self, value):
        self.well_upper_BHP.append(value)


def AssembleDefaultSimulationData(path_to_ctrl):
    #read all necessary data from a control file
    #yet to be implemented...

    sim_defaults = ctrl_sto()

    '''
    DEBUG VALUES FROM HERE ONWARDS
    '''
    # these are variables which should be provided in a control structure (they do not change during sim)
    well_names = ['Well_c', 'Well_n', 'Well_ne', 'Well_e', 'Well_se', 'Well_s', 'Well_sw', 'Well_w', 'Well_nw']
    well_l_bhp = [35, 35, 35, 35, 35, 35, 35, 35, 35]
    well_u_bhp = [120, 120, 120, 120, 120, 120, 120, 120, 120]

    
    for entry in well_names:
       sim_defaults.add_well_name(entry)
    for entry in well_l_bhp:
       sim_defaults.add_well_lower_BHP_val(entry)
    for entry in well_u_bhp:
       sim_defaults.add_well_upper_BHP_val(entry)

    sim_defaults.set_path(r'E:\Programming\IF_PPlant_Storage\testdata')
    sim_defaults.set_simulator('e300')
    sim_defaults.set_simulation_title('SYNTH_ANTI_CAES_IF')
    sim_defaults.set_restart_id(0)
    sim_defaults.set_simulator_path(r'C:\ecl\2017.2\bin\pc_x86_64')
    sim_defaults.set_ecl_log_flag(True)

    '''
    END OF DEBUG VALUES
    '''

    return sim_defaults


def entry_node(target_flow, tstep, tstepsize, op_mode):

    #this is the entry point for the geostorage coupling
    #if tstep == 1:
    sim_defaults = AssembleDefaultSimulationData('')

    if(sim_defaults.simulator == 'ECLIPSE' or sim_defaults.simulator == 'e300'):
        flowrate, pressure = RunECLIPSE(sim_defaults, target_flow, tstepsize, tstep, op_mode)
    elif sim_defaults.simulator == 'proxy':
        pass
        #implement later
    else:
        print('ERROR: simulator flag not understood. Is: ', sim_defaults.simulator)
    
    return (flowrate, pressure)


def RunECLIPSE(ctrl, target_flowrate, tstepsize, tstep, current_mode):

    #wrapper for the eclipse coupling
    print('######################################################################')
    print('Running storage simulation for timestep:\t', '%.0f'%tstep)
    print('timestep size:\t\t\t\t\t', tstepsize, '\t\ts')
    print('target storage flowrate:\t\t\t', '%.6f'%target_flowrate, '\tsm3/s')
    print('operational mode:\t\t\t\t', current_mode)
    print('----------------------------------------------------------------------')
    
    #('%.2f'%a) 
    
    # assembling current ecl data file
    reworkECLData(ctrl, tstep, tstepsize, target_flowrate, current_mode)
    #print('tstep: ', tstep, ' stepsize: ', tstepsize)
    #print('current operational mode: ', current_mode, ' target flowrate: ', target_flowrate)
    # executing eclipse
    ExecuteECLIPSE(ctrl, tstep)
    # reading results
    ecl_results = GetECLResults(ctrl, tstep, current_mode)
    # returns achieved flowrate, pressure, ?

    print('Timestep ', tstep, ' completed.')
    print('Pressure actual:\t', '%.6f'%ecl_results[0], '\tbars')
    print('Flowrate actual:\t', '%.6f'%ecl_results[1], '\tsm3/s')
    print('######################################################################')
    return (ecl_results[1], ecl_results[0])



def writeFile(path, a_list):
    #writes a_list to a file
    with open (path, 'w') as f:
        for entry in a_list:
            f.write("%s" % entry)



def getFile(path):
    #opens file and returns it as list
    a_list = []
    with open(path) as f:
        a_list = list(f)
    f.closed
    
    return a_list



def rearrangeRSMDataArray(rsm_list):
    # function to re-order / get rid off line breaks etc. in input
    break_count = getStringCount(rsm_list, 'SUMMARY OF RUN')
    break_positions = getStringPositions(rsm_list, 'SUMMARY OF RUN')
    #break_positions = [i for i, s in enumerate(rsm_list) if 'SUMMARY OF RUN' in s]

    if break_count > 0:
        interval = break_positions[1] - break_positions[0]

    output = []

    for i in range(break_count):
        for j in range(interval):
            current_idx = i * interval + j
            if i == 0:
                temp = rsm_list[current_idx][:-3]
                output.append(rsm_list[i])
            if j > 1:
                #print('current_idx: ', current_idx)
                temp1 = str(output[j]).replace('\n', '')
                temp1 = temp1[:-3]
                temp2 = str(rsm_list[current_idx])
                temp2 = temp2[:-3]
                temp2 += '\n'
                temp = temp1 + temp2
                temp = temp.replace('\t\t', '\t')
                
                output[j] = temp
    #delete first two (empty) entries
    del output[0]
    del output[0]

    return output



def contractDataArray(input):
    # takes a list array as input
    # removes starting and intermediate whitespaces
    # returns a cleaner array 
    output = []
    rows = len(input)
    
    for i in range(rows):
        # for each row do remove the stuff and save as clean list
        a_row = input[i]
        a_row = a_row.replace('\t', ';')
        a_row = a_row.split(';')
       # if not a_row[0]: 
        del a_row[0]  #first entry is always blank in ecl rsm output
        a_new_row = []
        for j in a_row:
            temp_str = j.strip()
            if not str(temp_str):
                temp_str = 'n.a.'
            a_new_row.append(temp_str)
        # if itis the last row delete last entry of that row
        #print('i: ', i , ' total rows: ', rows)
        if (i + 1) == rows:
            del a_new_row[-1]
            
        output.append(a_new_row)
    
    #delete last entry in last row
    
    #go through whole list and delte double date entries
    date_positions = [i for i, s in enumerate(output[0]) if 'DATE' in s]
    #print('Positions of \'DATE\' in array: ', date_positions)
    max_count = len(date_positions)
    for e in reversed(date_positions):
        if max_count > 1:
            for i in range(len(output)):
                del output[i][e]
        max_count = max_count - 1
    
    return output



def searchSection(ecl_data_list, section):
    #function searches through list and returns position of a section
    #what if there is a whitespace behind keyword?
    pos = 0

    if section in ecl_data_list:
        pos = ecl_data_list.index(section)
    else:
        section += "\n"
        if section in ecl_data_list:
            pos = ecl_data_list.index(section)
    return pos



def getStringPositions(input, keyword):
    return [i for i, s in enumerate(input) if keyword in s]



def getStringCount(input, keyword):
    return sum ( 1 for s in input if keyword in s)



def reworkECLData(control, timestep, timestepsize, flowrate, op_mode):
        
    # open and read eclipse data file
    ecl_data_file = getFile(control.path + '\\' + control.simulation_title + '.DATA')
    
    #rearrange the entries in the saved list
    #first search for the restart keyword
    
    if timestep == 2:
        #look for EQUIL and RESTART keyword
        equil_pos = searchSection(ecl_data_file, 'EQUIL')
        if(equil_pos > 0):
            #delete equil and replace with restart
            #assemble new string for restart section
            ecl_data_file[equil_pos] = 'RESTART\n'
            ecl_data_file[equil_pos + 1] =  '\'' + control.simulation_title + '\' \t' 
            ecl_data_file[equil_pos + 1] += str(int(control.restart_id) + timestep - 1)  + ' /\n'
        else:
            restart_pos = searchSection(ecl_data_file, "RESTART")
            if restart_pos > 0:
                #assemble new string for restart section
                ecl_data_file[restart_pos + 1] =  '\'' + control.simulation_title + '\' \t' 
                ecl_data_file[restart_pos + 1] += str(int(control.restart_id) + timestep - 1)  + ' /\n'
    if timestep > 2:
        restart_pos = searchSection(ecl_data_file, "RESTART")
        if restart_pos > 0:
            #assemble new string for restart section
            ecl_data_file[restart_pos + 1] =  '\'' + control.simulation_title + '\' \t' 
            ecl_data_file[restart_pos + 1] += str(int(control.restart_id) + timestep - 1)  + ' /\n'
        
        
    #now rearrange the well schedule section
    schedule_pos = searchSection(ecl_data_file, "WCONINJE")
    if schedule_pos == 0:
        schedule_pos = searchSection(ecl_data_file, "WCONPROD")
        
    if schedule_pos > 0:
        # delete the old well schedule
        del ecl_data_file[schedule_pos:]
        # append new well schedule
        # first calculate rate applied for each well
        well_count = len(control.well_names)
        well_target = abs(flowrate / well_count)
        well_target_days = well_target * 60.0 * 60.0 *24.0

        #now construct new well schedule section
        ecl_data_file.append('\n')
        
        if op_mode == 'charging':
            ecl_data_file.append("WCONINJE\n")
            for idx in range(len(control.well_names)):
                line = '\'' + control.well_names[idx] + '\''
                line += '\t\'GAS\'\t\'OPEN\'\t\'RATE\'\t'
                line += str(well_target_days) + '\t'
                line += '1*\t' + str(control.well_upper_BHP[idx]) + '/\n'
                ecl_data_file.append(line)

        elif op_mode == 'discharging':
            ecl_data_file.append("WCONPROD\n")
            for idx in range(len(control.well_names)):
                line = '\'' + control.well_names[idx] + '\''
                line += '\t\'OPEN\'\t\'GRAT\'\t1*\t1*\t'
                line += str(well_target_days) + '\t'
                line += '1*\t1*\t' + str(control.well_lower_BHP[idx]) + '/\n'
                ecl_data_file.append(line)
        elif op_mode == 'shut-in':
            ecl_data_file.append("WCONINJE\n")
            for idx in range(len(control.well_names)):
                line = '\'' + control.well_names[idx] + '\''
                line += '\t\'GAS\'\t\'OPEN\'\t\'RATE\'\t'
                line += '0.0' + '\t'
                line += '1*\t' + str(control.well_upper_BHP[idx]) + '/\n'
                ecl_data_file.append(line)
        else:
            print('ERROR: operational mode not understood in timestep: ', timestep, ' is: ', op_mode)

        ecl_data_file.append('/')
        #finish schedule
        timestepsize_days = timestepsize / 60.0 / 60.0 / 24.0
        file_finish = ['\n', '\n', 'TSTEP\n', '1*' + str(timestepsize_days) + '\n', '/\n', '\n', '\n', 'END\n' ]
        ecl_data_file += file_finish

        #save to new file
        temp_path = control.path + '\\' + control.simulation_title + '.DATA'
        writeFile(temp_path, ecl_data_file)



def ExecuteECLIPSE(control, tstep):

    #import subprocess
    import os

    if os.name == 'nt':
        simulation_path = control.path + '\\' + control.simulation_title + '.DATA'
        if control.keep_ecl_logs == True:
            log_file_path = control.path + '\\' + 'log_' + control.simulation_title + '_' + str(tstep) + '.txt'
        else:
            log_file_path = 'NUL'
        temp = 'eclrun ' + control.simulator + ' ' + simulation_path + ' >' + log_file_path
        os.system(temp)
    #elif os.name == 'posix':
        #log_output_loc += ' &'
        #rc = subprocess.call(['eclrun', simulator_loc, simulation_title_loc, '>', log_output_loc])
    
    #print(rc)



def GetECLResults(control, timestep, current_op_mode):
    #function returns list with two entries:
    # data[pressure_actual, flowrate_actual]
    
    #first read the results file
    filename = control.path + '\\' + control.simulation_title + '.RSM'  
    results = getFile(filename)
    #sort the rsm data to a more uniform dataset
    reorderd_rsm_data = rearrangeRSMDataArray(results)
    #eleminate additional whitespaces, duplicate entries, etc.
    well_results = contractDataArray(reorderd_rsm_data)

    # check number of data entries in well_results:
    values = len(well_results) - 4
    if values > 1:
        print('Warning: possible loss of data, expecting one data line in RSM file only')
    
    #data structures to save the flowrates, pressures and names of all individual wells
    well_pressures = []
    well_flowrates_days = []
    well_flowrates = []
    well_mass_flowrates = [] #not in use yet
    well_names = []
    well_names_loc = []
    flowrate_actual = 0.0
    pressure_actual = 0.0
    
    # get well pressures
    pressure_keyword = 'WBHP'
    bhp_positions = getStringPositions(well_results[0], pressure_keyword)

    for i in bhp_positions:
        well_pressures.append(float(well_results[-1][i]))
        well_names.append(well_results[2][i])
    
    # now get well flowrates
    if current_op_mode == 'discharging':
        #negativ flowrates
        #get all positions of WGPR entries in well_results
        flow_keyword = 'WGPR'
        flow_positions = getStringPositions(well_results[0], flow_keyword)
        for i in flow_positions:
            well_flowrates_days.append(float(well_results[-1][i]) * -1.0) #ecl flowrates are always positive
            well_names_loc.append(well_results[2][i])
               
    elif current_op_mode == 'charging':
        #positive flowrates
        flow_keyword = 'WGIR'
        flow_positions = getStringPositions(well_results[0], flow_keyword)
        for i in flow_positions:
            well_flowrates_days.append(float(well_results[-1][i]))
            well_names_loc.append(well_results[2][i])

    elif current_op_mode == 'shut-in':
        #do nothing
        pass
    else: 
        print('Warning: operational mode not understood, assuming shut-in at timestep: ', timestep) 
     
    if ( current_op_mode == 'charging' or current_op_mode == 'discharging'):
        # go through well names list and compare strings.
        # rearrange if necessary to get correct match for pressures and flowrates
        correct_idx = []
        for i in range(len(well_names)):
            if well_names[i] == well_names_loc[i]:
                correct_idx.append(i)
            else:
                target_str = well_names[i]
                for j in range(len(well_names_loc)):
                    if well_names_loc[j] == target_str:
                        correct_idx.append(j)
        #sort entries in well_flowrates based on correct_idx
        well_flowrates_temp = well_flowrates_days
        for i in correct_idx:
            well_flowrates_days[i] = well_flowrates_temp[i]
    
        #calculate total flowrate 
        #maybe add timestep dependence (only needed if more than one per call)?
        #change unit of flowrates to sm3/s from sm3/d
        for i in range(len(well_flowrates_days)):
            well_flowrates.append(well_flowrates_days[i] / 60.0 / 60.0 / 24.0)

        flowrate_actual = sum(well_flowrates)

        #calculate average pressure
        pressure_actual = 0.0
        for i in range(len(well_pressures)):
            pressure_actual += well_pressures[i] * well_flowrates[i]
        pressure_actual = pressure_actual / flowrate_actual
    else:
        pressure_actual = sum(well_pressures) / float(len(well_pressures))

    # to-do: return mass flow rate instead of volumetric flow data
    return [pressure_actual, flowrate_actual]

'''debug'''
data = [0.0, 0.0]
data = entry_node(1.15741, 1, 3600, 'charging')
data = entry_node(0.0, 2, 3600, 'shut-in')
data = entry_node(-1.15741, 3, 3600, 'discharging')