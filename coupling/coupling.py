#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

__author__ = "witte, wtp"

"""

import pandas as pd
import numpy as np
import powerplant as pp
import geostorage as gs
import json
import datetime
import os

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
    path = (r'D:\Simulations\if_testcase\testcase.main_ctrl.json')
    cd = coupling_data(path=path)

    # create instances for power plant and storage
    geostorage = gs.geo_sto(cd)

    #min_well_depth = min(geostorage.well_depths)
    min_well_depth = 700

    powerplant = pp.model(cd, min_well_depth)


    input_ts = read_series(cd.working_dir + cd.input_timeseries_path)
    output_ts = pd.DataFrame(columns=['pressure', 'massflow',
                                      'massflow_actual', 'power_actual',
                                      'success'])

    '''debug values from here onwards'''
    #data = [0.0, 0.0]
    #data = geostorage.CallStorageSimulation(1.15741, 1, cd, 'charging')
    #data = geostorage.CallStorageSimulation(0.0, 2, cd, 'shut-in')
    #data = geostorage.CallStorageSimulation(-1.15741, 3, cd, 'discharging')
    '''end of debug values'''

    p0 = 0.0 #old pressure (from last time step / iter)
    # get initial pressure before the time loop
    p0, dummy_flow = geostorage.CallStorageSimulation(0.0, 0, 0, cd, 'init')

    last_time = cd.t_start
    for t_step in range(cd.t_steps_total):
        current_time = datetime.timedelta(seconds=t_step * cd.t_step_length) + cd.t_start

        try:
            target_power = input_ts.loc[current_time].power / 100 * 1e6
            last_time = current_time
        except KeyError:
            target_power = input_ts.loc[last_time].power / 100 * 1e6

        # calculate pressure, mass flow and power
        p, m, m_corr, power, success = calc_timestep(
                powerplant, geostorage, target_power, p0, cd, t_step)
        
        # save last pressure (p1) for next time step as p0
        p0 = p
        #deleting old files
        geostorage.deleteFFile(t_step)

        # write pressure, mass flow and power to .csv
        output_ts.loc[current_time] = np.array([p, m, m_corr, power, success])

        if t_step % cd.save_nth_t_step == 0:
            output_ts.to_csv(cd.working_dir + cd.output_timeseries_path)


def calc_timestep(powerplant, geostorage, power, p0, md, tstep):
    """
    calculates one timestep of coupled power plant - storage simulation

    :param powerplant: powerplant model
    :type powerplant: powerplant.model object
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
    tstep_accepted = False
    storage_mode = ''

    if power == 0:
        m = 0
        storage_mode = 'shut-in'
    elif power < 0:
        storage_mode = 'discharging'
        m = powerplant.get_mass_flow(power, p0, storage_mode)
    else:
        storage_mode = 'charging'
        m = powerplant.get_mass_flow(power, p0, storage_mode)

    if m == 0:
        storage_mode = "shut-in"


    #moved inner iteration into timestep function,
    #iterate until timestep is accepted
    p0_temp = p0

    for iter_step in range(md.max_iter): #do time-specific iterations
       

        if tstep_accepted:
            print('Message: Timestep accepted after iteration ', iter_step - 1)
            break
        
        print('Current iteration: ', iter_step, 'Mode is: ', storage_mode)

        #get pressure for the given target rate and the actually achieved flow rate from storage simulation
        p1, m_corr = geostorage.CallStorageSimulation(m, tstep, iter_step, md, storage_mode )
        

        if storage_mode == 'charging' or storage_mode == 'discharging':
            # pressure check
            if abs((p0_temp - p1) / p1) > md.pressure_diff_rel or abs((p0_temp- p1)) > md.pressure_diff_abs:
                m = powerplant.get_mass_flow(power, p1, storage_mode)
                print('Adjusting mass flow rate! ', '\t', p0_temp, p1, m, m_corr)
            # if pressure check is successful, mass flow check:
            # check for difference due to pressure limitations
            elif abs(m_corr) <= 1E-5:
                power = 0
                tstep_accepted = True
                print('Adjusting power to zero! ', '\t', p0_temp, p1, m, m_corr)

            elif abs((m - m_corr) / m_corr) > md.flow_diff_rel:
                power = powerplant.get_power(m_corr, p1, storage_mode)
                tstep_accepted = True
                print('Adjusting power to ', power, ' ! ', '\t', p0_temp, p1, m, m_corr)

            else:
                #return p1, m_corr, power
                tstep_accepted = True
                m = m_corr
        
        elif storage_mode == "shut-in":
            print('Force accepting timestep b/c storage shut-in')
            tstep_accepted = True
        else:
            print('Problem: Storage mode not understood')
            tstep_accepted = True
        #saving old pressure
        p0_temp = p1

    if not tstep_accepted:
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
        self.path = path

        with open(path) as f:
            self.__dict__.update(json.load(f))


        self.coupled_simulation()

    def coupled_simulation(self):
        """
        Function to set all required default data, e.g. well names, paths, ...

        :returns: no return value
        """

        str_tmp = self.path.strip('.main_ctrl.json')
        self.scenario = ""
        self.working_dir = ""

        i = 0
        key = ""
        if os.name == 'nt':
            key = "\\"
        elif os.name == 'posix':
            key = "/"
        else:
            print('Error: OS not supported')


        for c in str_tmp[::-1]:
            if c == key:
                self.working_dir = str_tmp[:-i]
                break
            self.scenario += c
            i += 1

        self.scenario = self.scenario[::-1]
        self.debug = bool(self.debug)
        date_format = '%Y-%m-%d %H:%M:%S'
        self.t_start = datetime.datetime.strptime(self.t_start, date_format)

        print('Reading inputile \"' + self.scenario + '.main_ctrl.json\" '
              'in working directory \"' + self.working_dir + '\"')

'''        if self.debug:
            print('DEBUG-OUTPUT for main control data')
            print('Time series path:\t' + self.input_timeseries_path)
            print('Start time:\t' + str(self.t_start))
            print('Time step length:\t' + str(self.t_step_length))
            print('Number of time steps:\t' + str(self.t_steps_total))
            print('Iteration limits:\t' + str(self.min_iter) + '\t' +
                  str(self.max_iter))
            print('Pressure convergence criteria:\t' +
                  str(self.pressure_diff_abs) +
                  ' bars\t' + str(self.pressure_diff_rel * 100) + ' %')
            print('END of DEBUG-OUTPUT for main control data')'''


__main__()
