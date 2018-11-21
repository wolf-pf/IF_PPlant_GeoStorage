#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@author: wtp
"""

class geo_sto:
    
    '''
    class to include geologic storage simulations

    '''
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


    def InitializeStorageDefaults(self, path_to_ctrl):
        '''
        Function to set all required default data, e.g. well names, paths, ...

        :param path_to_ctrl: path to input file containing control data
        :param type: str
        :returns: no return value
        '''

        # open and read control file
        sto_ctrl_list = getFile(path_to_ctrl)
        
        for i in range(len(sto_ctrl_list)):
            sto_ctrl_list[i] = sto_ctrl_list[i].strip()
        #print(sto_ctrl_list)
        sto_ctrl_list = list(filter(None, sto_ctrl_list))
    
        #search for simulation data
        #get beginning and end 
        if searchSection(sto_ctrl_list, '<storage_simulator_data>') >= 0:
            pos_ident = searchSection(sto_ctrl_list, '<identifier>')
            end_pos_ident = searchSection(sto_ctrl_list, '</identifier>')
            if pos_ident < end_pos_ident and pos_ident >= 0:
                self.set_simulator(sto_ctrl_list[pos_ident + 1])
    
            pos_simulator_path = searchSection(sto_ctrl_list, '<simulator_path>')
            end_pos_simulator_path = searchSection(sto_ctrl_list, '</simulator_path>')
            if pos_simulator_path < end_pos_simulator_path and pos_simulator_path >= 0:
                self.set_simulator_path(sto_ctrl_list[pos_simulator_path + 1])
            
            pos_sim_title = searchSection(sto_ctrl_list, '<simulation_title>')
            end_pos_sim_title = searchSection(sto_ctrl_list, '</simulation_title>')
            if pos_sim_title < end_pos_sim_title and pos_sim_title >= 0:
                self.set_simulator_path(sto_ctrl_list[pos_sim_title + 1])
    
            pos_path = searchSection(sto_ctrl_list, '<simulation_path>')
            end_pos_path = searchSection(sto_ctrl_list, '</simulation_path>')
            if pos_path < end_pos_path and pos_path >= 0:
                self.set_path(sto_ctrl_list[pos_path + 1])
    
            pos_res_id = searchSection(sto_ctrl_list, '<restart_id>')
            end_pos_res_id = searchSection(sto_ctrl_list, '</restart_id>')
            if pos_res_id < end_pos_res_id and pos_res_id >= 0:
                self.set_restart_id(sto_ctrl_list[pos_res_id + 1])
            
            pos_ecl_log_flag = searchSection(sto_ctrl_list, '<ecl_log_file>')
            end_pos_ecl_log_flag = searchSection(sto_ctrl_list, '</ecl_log_file>')
    
            if pos_ecl_log_flag < end_pos_ecl_log_flag and pos_ecl_log_flag >= 0:
                flag = sto_ctrl_list[pos_ecl_log_flag + 1]
                if bool(flag) == True or bool(flag) == False:
                    self.set_ecl_log_flag(bool(flag))
                else: 
                    print('Input for ecl log file not understood. Setting flag to True')
                    self.set_ecl_log_flag(True)
        
        if searchSection(sto_ctrl_list, '<well_data>') >= 0:
            pos_well_count = searchSection(sto_ctrl_list, '<count>')
            end_pos_well_count = searchSection(sto_ctrl_list, '</count>')

            if pos_well_count < end_pos_well_count and pos_well_count >= 0:
                well_count_loc = sto_ctrl_list[pos_well_count + 1]

                pos_well_data = searchSection(sto_ctrl_list, '<data>')
                end_pos_well_data = searchSection(sto_ctrl_list, '</data>')
                if pos_well_data < end_pos_well_data and pos_well_data >= 0:
                    for i in range(int(well_count_loc)):
                        temp_list = (sto_ctrl_list[pos_well_data + 1 + i]).split()
                        self.add_well_name(temp_list[0])
                        self.add_well_lower_BHP_val(float(temp_list[1]))
                        self.add_well_upper_BHP_val(float(temp_list[2]))
                        #print(' Well name: ', temp_list[0], '\tBHP range: ', temp_list[1], '-', temp_list[2], ' bars')


    def CallStorageSimulation(self, target_flow, tstep, tstepsize, op_mode):
        '''
        Entry point for geo-storage simulation, handles all data transfer, executes simulator
        and provides simulation results to power plant simulator

        :param target_flow: target storage flow rate in sm3/d; 15.5556°C, 1 atm
        :param type: float
        :param tstep: current timestep
        :param type: int
        :param tstepsize: length of current timestep
        :param type: float
        :param op_mode: current operational mode, either 'charging', 'discharging' or 'shut-in' 
        :param type: str
        :returns: returns tuple of new pressure at the well (in reservoir) and actual (achieved) storage flow rate
        '''
        #this is the entry point for the geostorage coupling
    
        if(self.simulator == 'ECLIPSE' or self.simulator == 'e300'):
            flowrate, pressure = self.RunECLIPSE(target_flow, tstep, tstepsize, op_mode)
        elif self.simulator == 'proxy':
            pass
            #implement later
        else:
            print('ERROR: simulator flag not understood. Is: ', self.simulator)
        
        return (flowrate, pressure)
    
    
    def RunECLIPSE(self, target_flowrate, tstep, tstepsize, current_mode):
        '''
        Function acting as a wrapper for using eclipse (SChlumberger) as a storage simulator

        :param target_flowrate: target storage flow rate in sm3/d; 15.5556°C, 1 atm
        :param type: float
        :param tstep: current timestep
        :param type: int
        :param tstepsize: length of current timestep
        :param type: float
        :param op_mode: current operational mode, either 'charging', 'discharging' or 'shut-in' 
        :param type: str
        :returns: returns tuple of new pressure at the well (in reservoir) and actual (achieved) storage flow rate
        '''

        print('######################################################################')
        print('Running storage simulation for timestep:\t', '%.0f'%tstep)
        print('timestep size:\t\t\t\t\t', tstepsize, '\t\ts')
        print('target storage flowrate:\t\t\t', '%.6f'%target_flowrate, '\tsm3/s')
        print('operational mode:\t\t\t\t', current_mode)
        print('----------------------------------------------------------------------')
        
        # assembling current ecl data file
        self.reworkECLData(tstep, tstepsize, target_flowrate, current_mode)
        # executing eclipse
        self.ExecuteECLIPSE(tstep)
        # reading results
        ecl_results = self.GetECLResults(tstep, current_mode)
        # returns achieved flowrate, pressure, ?
    
        print('Timestep ', tstep, ' completed.')
        print('Pressure actual:\t', '%.6f'%ecl_results[0], '\tbars')
        print('Flowrate actual:\t', '%.6f'%ecl_results[1], '\tsm3/s')
        print('######################################################################')
        return (ecl_results[1], ecl_results[0])
    

    def rearrangeRSMDataArray(self, rsm_list):
        '''
        Function to sort through Eclipse's RSM file and obtain well data from last timestep

        :param rsm_list: list containing the RSM file
        :param type: str
        :returns: returns a clearer version of the input list (type: list of strings)
        '''
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
    
    
    def reworkECLData(self, timestep, timestepsize, flowrate, op_mode):
        '''
        function to change settings in the eclipse input file required for the storage simulation

        :param timestep: current timestep of simulation
        :param type: int
        :param timestepsize: length of current timestep
        :param type: float
        :param flowrate: current target storage flow rate from power plant simulation
        :param type: float
        :param op_mode: current operational mode, either 'charging', 'discharging' or 'shut-in' 
        :param type: str
        :returns: no return value
        '''
        # open and read eclipse data file
        ecl_data_file = getFile(self.path + '\\' + self.simulation_title + '.DATA')
        
        #rearrange the entries in the saved list
        if timestep == 2:
            #look for EQUIL and RESTART keyword
            equil_pos = searchSection(ecl_data_file, 'EQUIL')
            if(equil_pos > 0):
                #delete equil and replace with restart
                #assemble new string for restart section
                ecl_data_file[equil_pos] = 'RESTART\n'
                ecl_data_file[equil_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                ecl_data_file[equil_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
            else:
                restart_pos = searchSection(ecl_data_file, "RESTART")
                if restart_pos > 0:
                    #assemble new string for restart section
                    ecl_data_file[restart_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                    ecl_data_file[restart_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
        if timestep > 2:
            restart_pos = searchSection(ecl_data_file, "RESTART")
            if restart_pos > 0:
                #assemble new string for restart section
                ecl_data_file[restart_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                ecl_data_file[restart_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
            
        #now rearrange the well schedule section
        schedule_pos = searchSection(ecl_data_file, "WCONINJE")
        if schedule_pos == 0:
            schedule_pos = searchSection(ecl_data_file, "WCONPROD")
            
        if schedule_pos > 0:
            # delete the old well schedule
            del ecl_data_file[schedule_pos:]
            # append new well schedule
            # first calculate rate applied for each well
            well_count = len(self.well_names)
            well_target = abs(flowrate / well_count)
            well_target_days = well_target * 60.0 * 60.0 *24.0
    
            #now construct new well schedule section
            ecl_data_file.append('\n')
            
            if op_mode == 'charging':
                ecl_data_file.append("WCONINJE\n")
                for idx in range(len(self.well_names)):
                    line = '\'' + self.well_names[idx] + '\''
                    line += '\t\'GAS\'\t\'OPEN\'\t\'RATE\'\t'
                    line += str(well_target_days) + '\t'
                    line += '1*\t' + str(self.well_upper_BHP[idx]) + '/\n'
                    ecl_data_file.append(line)
    
            elif op_mode == 'discharging':
                ecl_data_file.append("WCONPROD\n")
                for idx in range(len(self.well_names)):
                    line = '\'' + self.well_names[idx] + '\''
                    line += '\t\'OPEN\'\t\'GRAT\'\t1*\t1*\t'
                    line += str(well_target_days) + '\t'
                    line += '1*\t1*\t' + str(self.well_lower_BHP[idx]) + '/\n'
                    ecl_data_file.append(line)
            elif op_mode == 'shut-in':
                ecl_data_file.append("WCONINJE\n")
                for idx in range(len(self.well_names)):
                    line = '\'' + self.well_names[idx] + '\''
                    line += '\t\'GAS\'\t\'OPEN\'\t\'RATE\'\t'
                    line += '0.0' + '\t'
                    line += '1*\t' + str(self.well_upper_BHP[idx]) + '/\n'
                    ecl_data_file.append(line)
            else:
                print('ERROR: operational mode not understood in timestep: ', timestep, ' is: ', op_mode)
    
            ecl_data_file.append('/')
            #finish schedule
            timestepsize_days = timestepsize / 60.0 / 60.0 / 24.0
            file_finish = ['\n', '\n', 'TSTEP\n', '1*' + str(timestepsize_days) + '\n', '/\n', '\n', '\n', 'END\n' ]
            ecl_data_file += file_finish
    
            #save to new file
            temp_path = self.path + '\\' + self.simulation_title + '.DATA'
            writeFile(temp_path, ecl_data_file)
    
    
    
    def ExecuteECLIPSE(self, tstep):
        '''
        Function to call eclipse executable

        :param tstep: current timestep
        :param type: int
        :returns: no return value
        '''
        #import subprocess
        import os
    
        if os.name == 'nt':
            simulation_path = self.path + '\\' + self.simulation_title + '.DATA'
            if self.keep_ecl_logs == True:
                log_file_path = self.path + '\\' + 'log_' + self.simulation_title + '_' + str(tstep) + '.txt'
            else:
                log_file_path = 'NUL'
            temp = 'eclrun ' + self.simulator + ' ' + simulation_path + ' >' + log_file_path
            os.system(temp)
        #elif os.name == 'posix':
            #log_output_loc += ' &'
            #rc = subprocess.call(['eclrun', simulator_loc, simulation_title_loc, '>', log_output_loc])
        
        #print(rc)
    
    
    def GetECLResults(self, timestep, current_op_mode):
        '''
        Function to get the eclipse results from the *.RSM file and analyze the results

        :param timestep: current timestep
        :param type: int
        :param current_op_mode: operational mode, either 'charging', 'discharging' or 'shut-in'
        :param type: str
        :returns: returns a tuple of float values containing pressure and actual storage flow rate
        '''
        
        #first read the results file
        filename = self.path + '\\' + self.simulation_title + '.RSM'  
        results = getFile(filename)
        #sort the rsm data to a more uniform dataset
        reorderd_rsm_data = self.rearrangeRSMDataArray(results)
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
        
        # now get well flow rates
        if current_op_mode == 'discharging': #negative flow rates
            #get all positions of WGPR entries in well_results
            flow_keyword = 'WGPR'
            flow_positions = getStringPositions(well_results[0], flow_keyword)
            for i in flow_positions:
                well_flowrates_days.append(float(well_results[-1][i]) * -1.0) #ecl flowrates are always positive
                well_names_loc.append(well_results[2][i])
                   
        elif current_op_mode == 'charging': #positive flow rates
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
    
    

def writeFile(path, a_list):
    '''
    Short function to write a file based on a list of strings

    :param a_list: the input list
    :param type: string
    :returns: no return value
    '''
    with open (path, 'w') as f:
        for entry in a_list:
            f.write("%s" % entry)



def getFile(path):
    '''
    Short function to read a file and save as a list fo strings

    :param path: the path to the file 
    :param type: string
    :returns: a list of strings
    '''
    a_list = []
    with open(path) as f:
        a_list = list(f)
    f.closed
    
    return a_list


def contractDataArray(input):
    '''
    Short function to clean and contract a data array

    :param input: the input list to be cleaned
    :param type: string
    :returns: a list of strings
    '''
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
    
    #go through whole list and delete double date entries
    date_positions = [i for i, s in enumerate(output[0]) if 'DATE' in s]
    #print('Positions of \'DATE\' in array: ', date_positions)
    max_count = len(date_positions)
    for e in reversed(date_positions):
        if max_count > 1:
            for i in range(len(output)):
                del output[i][e]
        max_count = max_count - 1
    
    return output



def searchSection(data_list, section):
    '''
    function to search for a given string in a list of strings

    :param data_list: the input list to be searched
    :param type: list of strings
    :param section: the string which is searched for
    :param type: string
    :returns: int
    '''
    #what if there is a whitespace behind keyword?
    pos = -1
    if section in data_list:
        pos = data_list.index(section)
    else:
        section += "\n"
        if section in data_list:
            pos = data_list.index(section)
    return pos



def getStringPositions(input, keyword):
    '''
    function to get all positions of a string in a list

    :param input: the input list to be searched in
    :param type: list of strings
    :param keyword: the string to be searched for
    :param type: string
    :returns: list of int
    '''
    return [i for i, s in enumerate(input) if keyword in s]



def getStringCount(input, keyword):
    '''
    function to count the occurence of a string in a list

    :param input: the input list to be searched in
    :param type: list of strings
    :param keyword: the string to be searched for
    :param type: string
    :returns: int
    '''
    return sum ( 1 for s in input if keyword in s)

