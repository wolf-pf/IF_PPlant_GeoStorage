#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

__author__ = "wtp"

"""

import utilities as util

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
        #self.debug_flag = False

    #self.debug_flag = kwargs.get('DEBUG_FLAG')
    
    def set_ecl_log_flag(self, a_flag):
        self.keep_ecl_logs = a_flag

    #def set_debug_flag(self, a_flag):
    #    self.debug_flag = a_flag

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


    def InitializeStorageDefaults(self, path_to_ctrl, debug):
        '''
        Function to set all required default data, e.g. well names, paths, ...

        :param path_to_ctrl: path to input file containing control data
        :param type: str
        :returns: no return value
        '''

        # read and clean control file
        sto_ctrl_list = util.cleanControlFileList(util.getFile(path_to_ctrl))       
    
        #search for simulation data
        #get beginning and end 
        if util.getValuefromControlFileList(sto_ctrl_list, 'storage_simulator_data') is not 'KEY_NOT_FOUND':
            self.set_simulator(util.getValuefromControlFileList(sto_ctrl_list, 'identifier'))
            self.set_simulator_path(util.getValuefromControlFileList(sto_ctrl_list, 'simulator_path'))
            self.set_simulation_title(util.getValuefromControlFileList(sto_ctrl_list, 'simulation_title'))
            self.set_path(util.getValuefromControlFileList(sto_ctrl_list, 'simulation_path'))
            self.set_restart_id(util.getValuefromControlFileList(sto_ctrl_list, 'restart_id'))
            flag = util.getValuefromControlFileList(sto_ctrl_list, 'ecl_log_file')
            if bool(flag) == True or bool(flag) == False:
                self.set_ecl_log_flag(bool(flag))
            else: 
                print('Input for ecl log file not understood. Setting flag to True')
                self.set_ecl_log_flag(True)
        else:
            print('ERROR: No storage simulator data found in file:')
            print(path_to_ctrl)
        
        if util.getValuefromControlFileList(sto_ctrl_list, 'well_data') is not 'KEY_NOT_FOUND':
            well_count_loc = util.getValuefromControlFileList(sto_ctrl_list, 'count')
            if well_count_loc is not 'KEY_NOT_FOUND':
                well_count_loc = int(well_count_loc)
                
                if util.getValuefromControlFileList(sto_ctrl_list, 'data') is not 'KEY_NOT_FOUND':
                    start_idx = util.getIdxfromControlFileList(sto_ctrl_list, 'data')
                    for i in range(well_count_loc):
                        temp_list = (sto_ctrl_list[start_idx + 1 + i]).split()
                        self.add_well_name(temp_list[0])
                        self.add_well_lower_BHP_val(float(temp_list[1]))
                        self.add_well_upper_BHP_val(float(temp_list[2]))
                        #print(' Well name: ', temp_list[0], '\tBHP range: ', temp_list[1], '-', temp_list[2], ' bars')
        else:
            print('ERROR: No well data found in file:')
            print(path_to_ctrl)

        if  debug == True:
            print('DEBUG-OUTPUT for storage control data')
            print('Selected simulator:\t', self.simulator)
            print('Simulator path:\t', self.simulator_path)
            print('Simulation title:\t', self.simulation_title)
            print('Simulation path:\t', self.path)
            print('Restart ID:\t', self.restart_id)
            print('Safe ECL log files:\t', self.keep_ecl_logs)
            print('Number of wells:\t', str(len(self.well_names)))
            for i in range(len(self.well_names)):
                print('\t', self.well_names[i], ', ', self.well_lower_BHP[i], 'bars, ', self.well_upper_BHP[i], 'bars')
            print('END of DEBUG-OUTPUT for storage control data')



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
        
        return pressure, flowrate
    
    
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
        if not current_mode == 'init':
            print('Running storage simulation for timestep:\t', '%.0f'%tstep)
            print('timestep size:\t\t\t\t\t', tstepsize, '\t\ts')
            print('target storage flowrate:\t\t\t', '%.6f'%target_flowrate, '\tsm3/s')
            print('operational mode:\t\t\t\t', current_mode)
            print('----------------------------------------------------------------------')
        else:
            print('Running storage simulation to obtain initial pressure')
        
        # assembling current ecl data file
        self.reworkECLData(tstep, tstepsize, target_flowrate, current_mode)
        # executing eclipse
        self.ExecuteECLIPSE(tstep, current_mode)
        # reading results
        ecl_results = self.GetECLResults(tstep, current_mode)
        # returns achieved flowrate, pressure, ?
    
        if not current_mode == 'init':
            print('Timestep ', tstep, ' completed.')
            print('Pressure actual:\t', '%.6f'%ecl_results[0], '\tbars')
            print('Flowrate actual:\t', '%.6f'%ecl_results[1], '\tsm3/s')
        else:
            print('Running storage simulation to obtain initial pressure')
            print('Initial pressure is: \t', '%.6f'%ecl_results[0], '\tbars')
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
        break_count = util.getStringCount(rsm_list, 'SUMMARY OF RUN')
        break_positions = util.getStringPositions(rsm_list, 'SUMMARY OF RUN')
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
        ecl_data_file = util.getFile(self.path + '\\' + self.simulation_title + '.DATA')
        
        #rearrange the entries in the saved list
        if timestep == 2:
            #look for EQUIL and RESTART keyword
            equil_pos = util.searchSection(ecl_data_file, 'EQUIL')
            if(equil_pos > 0):
                #delete equil and replace with restart
                #assemble new string for restart section
                ecl_data_file[equil_pos] = 'RESTART\n'
                ecl_data_file[equil_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                ecl_data_file[equil_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
            else:
                restart_pos = util.searchSection(ecl_data_file, "RESTART")
                if restart_pos > 0:
                    #assemble new string for restart section
                    ecl_data_file[restart_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                    ecl_data_file[restart_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
        if timestep > 2:
            restart_pos = util.searchSection(ecl_data_file, "RESTART")
            if restart_pos > 0:
                #assemble new string for restart section
                ecl_data_file[restart_pos + 1] =  '\'' + self.simulation_title + '\' \t' 
                ecl_data_file[restart_pos + 1] += str(int(self.restart_id) + timestep - 1)  + ' /\n'
            
        #now rearrange the well schedule section
        schedule_pos = util.searchSection(ecl_data_file, "WCONINJE")
        if schedule_pos == 0:
            schedule_pos = util.searchSection(ecl_data_file, "WCONPROD")
            
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

            elif op_mode == 'shut-in' or op_mode == 'init':
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
            if not op_mode == 'init':
                temp_path = self.path + '\\' + self.simulation_title + '.DATA'
            else:
                temp_path = self.path + '\\' + self.simulation_title + '_init.DATA'
            util.writeFile(temp_path, ecl_data_file)
    
    
    
    def ExecuteECLIPSE(self, tstep, op_mode):
        '''
        Function to call eclipse executable

        :param tstep: current timestep
        :param type: int
        :param op_mode: operational mode of storage simulation
        :param type: str
        :returns: no return value
        '''
        #import subprocess
        import os
    
        if os.name == 'nt':
            simulation_path = ''
            if not op_mode == 'init':
                simulation_path = self.path + '\\' + self.simulation_title + '.DATA'
            else:
                simulation_path = self.path + '\\' + self.simulation_title + '_init.DATA'

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
        if not current_op_mode == 'init':
            filename = self.path + '\\' + self.simulation_title + '.RSM'
        else:
            filename = self.path + '\\' + self.simulation_title + '_init.RSM'
        results = util.getFile(filename)
        #sort the rsm data to a more uniform dataset
        reorderd_rsm_data = self.rearrangeRSMDataArray(results)
        #eleminate additional whitespaces, duplicate entries, etc.
        well_results = util.contractDataArray(reorderd_rsm_data)
    
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
        bhp_positions = util.getStringPositions(well_results[0], pressure_keyword)
    
        for i in bhp_positions:
            well_pressures.append(float(well_results[-1][i]))
            well_names.append(well_results[2][i])
            if well_pressures[-1] == 0.0:
                print('Problem: well pressure for well ', well_names[-1], ' is zero. Setting to corresponding BHP limit' )
                bhp_limits_well = self.getWellBHPLimits(well_names[-1])
                if current_op_mode == 'discharging':
                    well_pressures[-1] = bhp_limits_well[0]
                elif current_op_mode == 'charging' or current_op_mode == 'shut-in':
                    well_pressures[-1] = bhp_limits_well[1]
                else:
                    print('Problem: could not determine operational mode, assuming injection')
                    well_pressures[-1] = bhp_limits_well[1]
                
        
        # now get well flow rates
        if current_op_mode == 'discharging': #negative flow rates
            #get all positions of WGPR entries in well_results
            flow_keyword = 'WGPR'
            flow_positions = util.getStringPositions(well_results[0], flow_keyword)
            for i in flow_positions:
                well_flowrates_days.append(float(well_results[-1][i]) * -1.0) #ecl flowrates are always positive
                well_names_loc.append(well_results[2][i])
                   
        elif current_op_mode == 'charging': #positive flow rates
            flow_keyword = 'WGIR'
            flow_positions = util.getStringPositions(well_results[0], flow_keyword)
            for i in flow_positions:
                well_flowrates_days.append(float(well_results[-1][i]))
                well_names_loc.append(well_results[2][i])
    
        elif current_op_mode == 'shut-in' or current_op_mode == 'init':
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
    
    def getWellBHPLimits(self, well_name):
        '''
        function to obtain pressure limits for a given well
        :param well_name: well identifier used to search well list
        :param type: string
        :returns: tuple of float, lower and upper BHP limit
        '''
        for i in range(len(self.well_names)):
            if self.well_names[i] == well_name:
                return [self.well_lower_BHP[i], self.well_upper_BHP[i]]

        return [0.0, 0.0]




