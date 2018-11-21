#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

__author__ = "witte, wtp"

"""

import pandas as pd
import numpy as np
import power_plant as pp
#import storage as sto
import geostorage as gs


def __main__(md):
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

    # create instances for power plant and storage
    power_plant = pp.model(md.pp)
    #storage = sto.model(md.sto)
    geostorage = gs.geo_sto() 

    time_series = read_series(md.ts_path)

    #get input control file for storage simulation
    '''debug values from here onwards'''
    #if tstep == 1:
    geostorage.InitializeStorageDefaults(r'.\testdata\testcase.storage_ctrl')
    #data = [0.0, 0.0]
    #data = CallStorageSimulation(1.15741, 1, 3600, 'charging')
    #data = CallStorageSimulation(0.0, 2, 3600, 'shut-in')
    #data = CallStorageSimulation(-1.15741, 3, 3600, 'discharging')
    '''end of debug values'''

    for tstep in time_series.index: #what is this loop doing? Is this the main time stepping loop?
        p0 = 0.0
        
        # get initial pressure (equals the pressure of the last timestep after the first iteration
        if tstep == 0: #do we start at 0 or 1?
            p0 = geostorage.CallStorageSimulation( 0.0, 0, md.stepwidth, 'init')
        
        # calculate pressure, mass flow and power
        p, m, power = calc_timestep(power_plant, geostorage,
                                        time_series.power[tstep ], p0, md, )
        #save last pressure (p1) for next time step as p0
        p0 = p
        
        # write pressure, mass flow and power to .csv
        write_timestep(md.out_path, np.array([tstep  + sub_step * md.stepwidth,
                                                  p, m, power]))


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

    for sub_step in range(md.num_timesteps): #are these sub-timesteps or time-specific iterations?
        if tstep_acctpted == True:
            break
        
        #get pressure for the given target rate and the actually achieved flow rate from storage simulation
        p1, m_corr = geostorage.CallStorageSimulation(m, tstep ,md.stepwidth, storage_mode )

        if storage_mode == 'charge' or storage_mode == 'discharge':
            # pressure check
            if ((abs((p0 - p1) / p1) > md.gap_rel or md.iter_count < md.min_iter) and md.iter_count < md.max_iter):
                md.iter_count += 1
                #calc_timestep(power_plant, storage, power, p1)
            # if pressure check is successful, mass flow check:
            # check for difference due to pressure limitations
            elif abs((m - m_corr) / m_corr) > md.gap_rel:
                power = power_plant.get_power(p1, m_corr)
                #return p1, m_corr, power_corr
            # return values
            else:
                #return p1, m_corr, power
                tstep_accepted = True
                m = m_corr
        else:
            tstep_accepted = True
        # return pressure if mass flow is zero
        #else:
        #    return p1, m, power
    if tstep_accepted != True:
        print('Problem: Results in timestep ', tstep, 'did not converge, accepting last iteration result.')
        
    return p1, m, power
    


def write_timestep(path, output):
    """
    writes the data of the timestep to .csv-file at path

    :param path: path to output file
    :type path: str
    :param output: output data for timestep
    :type output: list
    :returns: no return value
    """
    do_something


def read_series(path):
    """
    reads the input time series

    :param path: path to input time series
    :type path: str
    :returns: ts (*pandas.DataFrame*) - dataframe containing the time series
    """
    ts = pd.read_csv(path, index_col=0, delimiter=';', decimal='.')
    ts = ts.reindex(pd.to_datetime(ts.index))
    ts['power'] = ts['input_power'] - ts['output_power']
    return ts


class model_data:
    """
    creates a data container with the main model parameters

    :returns: no return value

    **allowed keywords** in kwargs:

    - ts (*str*) - path to input time series
    - out (*str*) - path for output time series
    - pp (*str*) - path to power plant data
    - sto (*str*) - path to storage data
    - min_iter (*int*) - minimum iterations on one timestep
    - max_iter (*int*) - maximum iterations on one timestep
    - gap_rel (*float*) - relative gap
    - gap_abs (*float*) - absoulte gap
    - stepwidth (*float*) - stepwidth as fractions of hours
    """

    def __init__(self, **kwargs):

        self.ts_path = kwargs.get('ts')
        self.out_path = kwargs.get('out')
        self.pp = kwargs.get('pp')
        self.sto = kwargs.get('sto')
        self.min_iter = kwargs.get('min_iter', 2)
        self.max_iter = kwargs.get('max_iter', 5)
        if self.min_iter >= self.max_iter:
            raise ValueError('Minimum number of iterations must be higher than'
                             ' than maximum number of iterations.')
        self.gap_rel = kwargs.get('gap_rel', 0.01)
        self.gap_abs = kwargs.get('gap_abs', 0.5)
        self.num_timesteps = kwargs.get('num_timesteps', 1)

        if not isinstance(self.num_timesteps, int):
            raise ValueError('Number of timesteps must be an integer value.')
        else:
            self.stepwidth = pd.Timedelta(3600 / self.num_timesteps, unit='s')

        self.iter_count = 0


md = model_data(ts='test.csv', out='results.csv',
                pp='power_plant/data.json', sto='storage/data.json')
__main__(md)
