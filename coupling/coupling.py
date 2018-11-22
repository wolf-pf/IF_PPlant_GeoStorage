#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

__author__ = "witte, wtp"

"""

import pandas as pd
import numpy as np
import power_plant as pp
import geostorage as gs
import utility as util
import json


def __main__():
    """
    main function to initialise the calculation

    - creates power plant and storage models
    - reads input timeseries
    - starts the loop
    - writes results to .csv-file

    :param md: object containing the basic model data
    :type md: model_data object
    :returns: no return value
    """
    #read main input file and set control variables, e.g. paths, identifiers, ...

    cd = coupling_data(path='testcase.main_ctrl.json')

    # create instances for power plant and storage
    power_plant = pp.model(md.pp)
    geostorage = gs.geo_sto()

    time_series = read_series(md.ts_path)
    ts_out = pd.Dataframe(columns=['pressure', 'massflow', 'massflow_actual',
                                   'power_actual', 'success'])

    #get input control file for storage simulation
    '''debug values from here onwards'''
    #if tstep == 1:
    geostorage.InitializeStorageDefaults(r'.\testdata\testcase.storage_ctrl', self.debug)
    #data = [0.0, 0.0]
    #data = CallStorageSimulation(1.15741, 1, 3600, 'charging')
    #data = CallStorageSimulation(0.0, 2, 3600, 'shut-in')
    #data = CallStorageSimulation(-1.15741, 3, 3600, 'discharging')
    '''end of debug values'''

    #for tstep in time_series.index: #what is this loop doing? Is this the main time stepping loop?
    for tstep_loc in range(self.tsteps_total): #what is this loop doing? Is this the main time stepping loop?
        current_time = tstep_loc * self.tstep_length + time_series.index[0] #anpassen!?

        target_power = time_series.loc[current_time].power / 100
        p0 = 0.0

        # get initial pressure (equals the pressure of the last timestep after the first iteration
        if tstep == 0: #do we start at 0 or 1?
            p0 = geostorage.CallStorageSimulation( 0.0, 0, self.tstep_length, 'init')

        # calculate pressure, mass flow and power
        p, m, m_corr, power, tstep_acctpted = calc_timestep(power_plant, geostorage,
                                        target_power, p0, md, tstep_loc )

        #save last pressure (p1) for next time step as p0
        p0 = p

        # write pressure, mass flow and power to .csv
        ts_out.loc[current_time] = np.array(
                [p, m, m_corr, power, tstep_acctpted])

        if tstep % 10 == 0:
            ts_out.to_csv(md.out_path)


    def InitializeCoupledSimulation(self, path_to_ctrl, debug):
        '''
        Function to set all required default data, e.g. well names, paths, ...
        :param path_to_ctrl: path to input file containing control data
        :param type: str
        :returns: no return value
        '''
        title_loc = ''
        dir_loc = ''

        #remove trailing slash if required
        if path_to_ctrl[-1] == '\\':
            del path_to_ctrl[-1]
        #get all slash positions for truncating name and path
        idxs = [i for i,key in enumerate(path_to_ctrl) if key=='\\']
        #append file ending if required
        if not path_to_ctrl[-10:] == '.main_ctrl':
            title_loc = path_to_ctrl[idxs[-1] + 1:]
            dir_loc = path_to_ctrl[:-len(title_loc)]
            path_to_ctrl += '\\.main_ctrl'
        else:
            pos_diff = len(path_to_ctrl) - idxs[-1] + 1
            title_loc = path_to_ctrl[idxs[-1] + 1:-(pos_diff - 10)]
            dir_loc = path_to_ctrl[:idxs[-1] + 1]


        print('Reading inputile: ', title_loc)
        print('in working directory: ', dir_loc)

        # read and clean control file
        main_ctrl_list = util.cleanControlFileList(util.getFile(path_to_ctrl))

        if util.getValuefromControlFileList(main_ctrl_list, 'main_simulation_data') is not 'KEY_NOT_FOUND':
            self.set_path_time_series(util.getValuefromControlFileList(main_ctrl_list, 'time_series_path'))
            self.set_tstep_length(util.getValuefromControlFileList(main_ctrl_list, 'time_step_length'))
            self.set_tsteps_total(util.getValuefromControlFileList(main_ctrl_list, 'time_steps'))
            temp_str = util.getValuefromControlFileList(main_ctrl_list, 'coupling_iterations')
            temp_lst = temp_str.split()
            if len(temp_lst) == 2:
                self.set_min_iter(temp_lst[0])
                self.set_max_iter(temp_lst[1])
            else:
                print('Problem: Coupling iteration keyword not understood, assuming 10 iterations max')
                self.set_min_iter(1)
                self.set_max_iter(10)
            self.set_max_abs_pdiff(util.getValuefromControlFileList(main_ctrl_list, 'max_absolute_pressure_diff'))
            self.set_max_rel_pdiff(util.getValuefromControlFileList(main_ctrl_list, 'max_relative_pressure_diff'))
        else:
            print('Problem: Keyword \'main_simulation_data\' not found in input file.')

        if debug == True:
            print('DEBUG-OUTPUT for main control data')
            print('Time series path:\t', self.path_time_series)
            print('Time step length:\t', self.tstep_length)
            print('Number of time steps:\t', self.tsteps_total)
            print('iteration limits:\t', self.min_iter, '\t', self.max_iter)
            print('Pressure convergence criteria:\t', self.max_pdiff_abs, 'bars\t', max_pdiff_rel, '%')
            print('END of DEBUG-OUTPUT for main control data')




def calc_timestep(power_plant, geostorage, power, p0, md, tstep):
    """
    calculates one timestep of coupled power plant - storage simulation

    :param power_plant: power plant model
    :type power_plant: power_plant.model object
    :param storage: storage model
    :type storage: storage.model object
    :param power: scheduled power for timestep
    :type power: float
    :param p0: initual pressure at timestep
    :type p0: float
    :param md: object containing the basic model data
    :type md: model_data object
    :returns: - p1 (*float*) - interface pressure at the end of the timestep
              - m_corr (*float*) - mass flow for this timestep
              - power (*float*) - power plant's input/output power for this
                timestep
    """
    tstep_acctpted = False
    storage_mode = ''

    if power == 0:
        m = 0
        storage_mode = 'shut-in'
    elif power < 0:
        storage_mode = 'discharge'
        m = power_plant.get_mass_flow(power, p0, storage_mode)
    else:
        storage_mode = 'charge'
        m = power_plant.get_mass_flow(power, p0, storage_mode)

    #moved inner iteration into timestep function,
    #iterate until timestep is accepted

    for iter_step in range(self.max_iter): #do time-specific iterations

        if tstep_acctpted == True:
            break

        #get pressure for the given target rate and the actually achieved flow rate from storage simulation
        p1, m_corr = geostorage.CallStorageSimulation(m, tstep ,self.tstep_length, storage_mode )

        if storage_mode == 'charge' or storage_mode == 'discharge':
            # pressure check
            if abs((p0 - p1) / p1) > md.gap_rel or abs((p0 - p1)) > md.gap_abs:
                continue
            # if pressure check is successful, mass flow check:
            # check for difference due to pressure limitations
            elif abs((m - m_corr) / m_corr) > md.gap_rel:
                power = power_plant.get_power(p1, m_corr)
                tstep_acctpted = True
            else:
                #return p1, m_corr, power
                tstep_accepted = True
                m = m_corr
        else:
            tstep_accepted = True

    if tstep_accepted != True:
        print('Problem: Results in timestep ', tstep, 'did not converge, accepting last iteration result.')

    return p1, m, m_corr, power, tstep_accepted


def read_series(path):
    """
    reads the input time series

    :param path: path to input time series
    :type path: str
    :returns: ts (*pandas.DataFrame*) - dataframe containing the time series
    """
    ts = pd.read_csv(path, delimiter=',', decimal='.')
    ts = ts.set_index(pd.to_datetime(ts.index, unit='h'))
    ts['power'] = ts['input'] - ts['output']
    return ts


class coupling_data:
    """
    creates a data container with the main model parameters

    :returns: no return value
    """

    def __init__(self, path):

        # load data.json information into objects dictionary (= attributes of
        # the object)
        with open(path) as f:
            self.__dict__.update(json.load(f))


__main__()
