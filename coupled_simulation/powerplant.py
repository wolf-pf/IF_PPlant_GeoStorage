#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

@author: witte
"""

# %% imports

import pandas as pd
import numpy as np
import json
from scipy import interpolate
import logging
from tespy import nwkr, logger, con, hlp

logger.define_logging(
    log_path=True, log_version=True, screen_level=logging.WARNING
)
# %% power plant model class


class model:
    """
    Creates the model for the power plant. Parameters are loaded from
    coupling data object cd.

    Parameters
    ----------
    cd : coupling_data
        Generel data for the interface handling.

    min_well_depth : float
        Depth of the wells.

    num_wells : int
        Number of wells.

    p_max : float
        Maximum pressure limit.

    p_min : float
        Minimum pressure limit.

    Note
    ----
    The depth of the wells along with the number of wells determines the
    dynamic pressure loss in bore hole pipes connecting the power plant with
    the geological storage. The pressure limits are the pressure limits at the
    bottom of the bore holes. These inforamtion are provided in
    the geological storage model control file.

    Example
    -------
    >>> from coupled_simulation import cp, pp
    """

    def __init__(self, cd, min_well_depth, num_wells, p_max, p_min):

        # load data.json information into objects dictionary (= attributes of
        # the object)

        # well information
        self.min_well_depth = min_well_depth
        self.num_wells = num_wells

        # pressure limits
        self.p_max = p_max
        self.p_min = p_min

        path = (cd.working_dir + cd.powerplant_path + cd.scenario + '.powerplant_ctrl.json')
        self.wdir = cd.working_dir + cd.powerplant_path
        with open(path) as f:
            self.__dict__.update(json.load(f))

        # load tespy models with the network_reader module
        self.tespy_charge = nwkr.load_nwk(self.wdir + self.tespy_charge_path)
        self.tespy_charge.set_printoptions(print_level='none')
        self.tespy_discharge = nwkr.load_nwk(self.wdir + self.tespy_discharge_path)
        self.tespy_discharge.set_printoptions(print_level='none')

        self.power_plant_layout()

# this part needs a rework
#        if self.spline_create_lut:
#            # create lookup table from tespy model
#            self.spline_charge_path = self.wdir + self.spline_charge_path
#            self.spline_discharge_path = self.wdir + self.spline_discharge_path
#
#        else:
#            # use existing lut at specified path
#            self.spline_charge_path = self.wdir + self.spline_charge_path
#            self.spline_discharge_path = self.wdir + self.spline_discharge_path
#
#        # load splines from .csv data
#        self.spline_charge = self.load_lookup(self.spline_charge_path)
#        self.spline_discharge = self.load_lookup(self.spline_discharge_path)

    def power_plant_layout(self):
        """
        Power plant layout calculation to determine power plant design point using
        nominal power input/output and nominal pressure as inputs.
        """
        # charging
        self.tespy_charge.imp_busses[self.power_bus_charge].set_attr(P=self.power_nominal_charge)
        self.tespy_charge.imp_conns[self.pressure_conn_charge].set_attr(
                p=self.pressure_nominal_charge,
                m=con.ref(self.tespy_charge.imp_conns[self.massflow_conn_charge], 1 / self.num_wells, 0))
        self.tespy_charge.imp_conns[self.massflow_conn_charge].set_attr(m=np.nan)
        self.tespy_charge.imp_comps[self.pipe_charge].set_attr(L=self.min_well_depth)
        self.tespy_charge.solve('design')
        self.tespy_charge.save(self.wdir + '/charge_design', path_abs=True)
        self.m_nom_charge = self.tespy_charge.imp_conns[self.massflow_conn_charge].m.val_SI
        self.m_min_charge = self.m_nom_charge * self.massflow_min_rel
        self.m_max_charge = self.m_nom_charge * self.massflow_max_rel
        msg = 'Nominal mass flow for charging is ' + str(self.m_nom_charge) + ' at nominal power ' + str(self.power_nominal_charge) + ' and nominal pressure ' + str(self.pressure_nominal_charge) + '.'
        logging.info(msg)

        # discharging
        self.tespy_discharge.imp_busses[self.power_bus_discharge].set_attr(P=self.power_nominal_discharge)
        self.tespy_discharge.imp_conns[self.pressure_conn_discharge].set_attr(
                p=self.pressure_nominal_discharge,
                m=con.ref(self.tespy_discharge.imp_conns[self.massflow_conn_discharge], 1 / self.num_wells, 0))
        self.tespy_discharge.imp_conns[self.massflow_conn_discharge].set_attr(m=np.nan)
        self.tespy_charge.imp_comps[self.pipe_discharge].set_attr(L=self.min_well_depth)
        self.tespy_discharge.solve('design')
        self.tespy_discharge.save(self.wdir + '/discharge_design', path_abs=True)
        self.tespy_discharge.print_results()
        self.m_nom_discharge = self.tespy_discharge.imp_conns[self.massflow_conn_discharge].m.val_SI
        self.m_min_discharge = self.m_nom_discharge * self.massflow_min_rel
        self.m_max_discharge = self.m_nom_discharge * self.massflow_max_rel
        msg = 'Nominal mass flow for discharging is ' + str(self.m_nom_discharge) + ' at nominal power ' + str(self.power_nominal_discharge) + ' and nominal pressure ' + str(self.pressure_nominal_discharge) + '.'
        logging.info(msg)

    def get_mass_flow(self, power, pressure, mode):
        """
        Calculate the mass flow at given power input (charging) or
        power output (discharging) and pressure at bottom borehole pressure.

        Parameters
        ----------
        power : float
            Scheduled electrical power input/output of the power plant.

        pressure : float
            Bottom borehole pressure.

        mode : str
            Calculation mode: :code:`mode in ['charging', 'discharging']`.

        Returns
        -------
        mass_flow : float
            Air mass flow from/into the storage.

        power_actual : float
            Actual electrical power input/output of the power plant.
            Differs from scheduled power, if schedule can not be met.
        """
        if pressure + 1e-4 < self.p_min:
            msg = 'Pressure is below minimum pressure: min=' + str(self.p_min) + 'value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0
        if pressure - 1e-4 > self.p_max:
            msg = 'Pressure is above maximum pressure: max=' + str(self.p_max) + 'value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0

        if self.method == 'tespy':
            if mode == 'charging':
                # if power too small
                if abs(power) < abs(self.power_nominal_charge / 100):
                    return 0, 0

                design_path = self.wdir + '/charge_design'
                # set power of bus
                self.tespy_charge.imp_busses[self.power_bus_charge].set_attr(P=power)
                # set pressure at interface
                self.tespy_charge.imp_conns[self.pressure_conn_charge].set_attr(p=pressure)
                self.tespy_charge.imp_conns[self.massflow_conn_charge].set_attr(m=np.nan)

                try:
                    self.tespy_charge.solve(mode='offdesign', design_path=design_path, path_abs=True)
                    m = self.tespy_charge.imp_conns[self.massflow_conn_charge].m.val_SI

                    if self.tespy_charge.res[-1] > 1e-3:
                        msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                        logging.error(msg)
                        return 0, 0
                    elif m < self.m_min_charge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' below minimum mass flow.'
                        logging.error(msg)
                        return 0, 0
                    elif m > self.m_max_charge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' above maximum mass flow. Adjusting power to match maximum allowed mass flow.'
                        logging.warning(msg)
                        return self.m_max_charge, self.get_power(self.m_max_charge, pressure, mode)
                    else:
                        msg = 'Calculation successful for power=' + str(power) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                        logging.debug(msg)
                        return m, power

                except:
                    # except general errors in calculation
                    msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                    logging.error(msg)
                    return 0, 0

            elif mode == 'discharging':
                if abs(power) < abs(self.power_nominal_discharge / 100):
                    return 0, 0

                design_path = self.wdir + '/discharge_design'
                # set power of bus
                self.tespy_discharge.imp_busses[self.power_bus_discharge].set_attr(P=power)
                # set pressure at interface
                self.tespy_discharge.imp_conns[self.pressure_conn_discharge].set_attr(p=pressure)
                self.tespy_discharge.imp_conns[self.massflow_conn_discharge].set_attr(m=np.nan)

                try:
                    self.tespy_discharge.solve(mode='offdesign', design_path=design_path, path_abs=True)
                    m = self.tespy_discharge.imp_conns[self.massflow_conn_discharge].m.val_SI
                    self.tespy_discharge.print_results()

                    if self.tespy_discharge.res[-1] > 1e-3:
                        msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                        logging.error(msg)
                        return 0, 0
                    elif m < self.m_min_discharge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' below minimum mass flow.'
                        logging.error(msg)
                        return 0, 0
                    elif m > self.m_max_discharge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' above maximum mass flow. Adjusting power to match maximum allowed mass flow.'
                        logging.warning(msg)
                        return self.m_max_discharge, self.get_power(self.m_max_discharge, pressure, mode)
                    else:
                        msg = 'Calculation successful for power=' + str(power) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                        logging.debug(msg)
                        return m, power

                except:
                    # except general errors in calculation
                    msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                    logging.error(msg)
                    return 0, 0

            else:
                raise ValueError('Mode must be charging or discharging.')

#        elif self.method == 'spline':
#            if mode == 'charging':
#                func = self.spline_charge
#
#            elif mode == 'discharging':
#                func = self.spline_discharge
#                power = -power
#            else:
#                raise ValueError('Mode must be charging or discharging.')
#
#            mass_flow = hlp.newton(hlp.reverse_2d, hlp.reverse_2d_deriv,
#                               [func, pressure, power], 0)
#
#            if mass_flow == 0:
#                print('ERROR: Could not find a solution for input pair: '
#                      'power=' + str(power) + ' pressure=' + str(pressure))
#
#            return mass_flow

        else:
            raise ValueError('Method must be tespy or spline.')

    def get_power(self, mass_flow, pressure, mode):
        """
        Calculates the power at given mass flow and pressure in charging or discharging mode.

        Parameters
        ----------
        mass_flow : float
            Mass flow.

        pressure : float
            Bottom borehole pressure.

        mode : str
            Calculation mode: :code:`mode in ['charging', 'discharging']`.

        Returns
        -------
        power : float
            Actual electrical power input/output of the power plant.
        """
        if pressure + 1e-4 < self.p_min:
            msg = 'Pressure is below minimum pressure: min=' + str(self.p_min) + 'value=' + str(pressure) + '.'
            logging.error(msg)
            return 0
        if pressure - 1e-4 > self.p_max:
            msg = 'Pressure is above maximum pressure: max=' + str(self.p_max) + 'value=' + str(pressure) + '.'
            logging.error(msg)
            return 0

        if self.method == 'tespy':
            if mode == 'charging':
                if mass_flow < self.m_min_charge - 1e-4:
                    msg = 'Mass flow is below minimum mass flow, shutting down power plant.'
                    logging.error(msg)
                    return 0
                elif mass_flow > self.m_max_charge + 1e-4:
                    msg = 'Mass flow above maximum mass flow. Adjusting mass flow to maximum allowed mass flow.'
                    logging.warning(msg)
                    return self.m_max_charge, self.get_power(self.m_max_charge, pressure, mode)


                design_path = self.wdir + '/charge_design'
                # set power of bus
                self.tespy_charge.imp_busses[self.power_bus_charge].set_attr(P=np.nan)
                # set pressure at interface
                self.tespy_charge.imp_conns[self.pressure_conn_charge].set_attr(p=pressure)
                self.tespy_charge.imp_conns[self.massflow_conn_charge].set_attr(m=mass_flow)

                self.tespy_charge.solve(mode='offdesign', design_path=design_path, path_abs=True)

                return self.tespy_charge.imp_busses[self.power_bus_charge].P.val

            elif mode == 'discharging':
                if mass_flow < self.m_min_discharge - 1e-4:
                    msg = 'Mass flow is below minimum mass flow, shutting down power plant.'
                    logging.error(msg)
                    return 0
                elif mass_flow > self.m_max_discharge + 1e-4:
                    msg = 'Mass flow above maximum mass flow. Adjusting mass flow to maximum allowed mass flow.'
                    logging.warning(msg)
                    return self.m_max_discharge, self.get_power(self.m_max_discharge, pressure, mode)


                design_path = self.wdir + '/discharge_design'
                # set power of bus
                self.tespy_discharge.imp_busses[self.power_bus_discharge].set_attr(P=np.nan)
                # set pressure at interface
                self.tespy_discharge.imp_conns[self.pressure_conn_discharge].set_attr(p=pressure)
                self.tespy_discharge.imp_conns[self.massflow_conn_discharge].set_attr(m=mass_flow)

                self.tespy_discharge.solve(mode='offdesign', design_path=design_path, path_abs=True)

                return self.tespy_discharge.imp_busses[self.power_bus_discharge].P.val

            else:
                raise ValueError('Mode must be charge or discharge.')

#        elif self.method == 'spline':
#            if mode == 'charging':
#                val = self.spline_charge.ev(massflow, pressure)
#                if abs((self.get_mass_flow(val, pressure, mode) - massflow) /
#                       massflow) < 1e-5:
#                    return val
#                else:
#                    return 0
#
#            elif mode == 'discharging':
#                val = self.spline_discharge.ev(massflow, pressure)
#                if abs((self.get_mass_flow(val, pressure, mode) - massflow) /
#                       massflow) < 1e-5:
#                    return val
#                else:
#                    return 0
#
#            else:
#                raise ValueError('Mode must be charge or discharge.')
        else:
            raise ValueError('Method must be tespy or spline.')

    def load_lookup(self, path):
        """
        Creates a rectangular bivariate spline object from data given in path.

        Parameters
        ----------
        path : str
            Relative path to the lookup table.

        Returns
        -------
        func : scipy.interpolate.RectBivariateSpline
            Spline interpolation object.
         """
        df = pd.read_csv(path, index_col=0)

        y = df.as_matrix()  # power

        x1 = df.index.get_values()  # mass flow
        if x1[0] > x1[-1]:
            x1 = x1[::-1]
            y = y[::-1]
        x2 = np.array(list(map(float, list(df))))  # pressure
        if x2[0] > x2[-1]:
            x2 = x2[::-1]
            y = y[:, ::-1]

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func
