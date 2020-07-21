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
from tespy.networks import load_network
from tespy.tools import logger
from tespy.connections import ref
from tespy.tools.helpers import TESPyNetworkError

logger.define_logging(
    log_path=True, log_version=True, screen_level=logging.WARNING, file_level=logging.WARNING
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
    """

    def __init__(self, cd, min_well_depth, num_wells, p_max, p_min):

        # load data.json information into objects dictionary (= attributes of
        # the object)
        path = (cd.working_dir + cd.powerplant_path + cd.scenario + '.powerplant_ctrl.json')
        self.wdir = cd.working_dir + cd.powerplant_path
        self.sc = cd.scenario
        with open(path) as f:
            self.__dict__.update(json.load(f))

        if self.method == 'tespy' or self.create_lut is True:
            # well information
            self.min_well_depth = min_well_depth
            self.num_wells = num_wells

            # pressure limits
            self.p_max = p_max
            self.p_min = p_min

            self.load_tespy_model()

            if self.create_lut is True:
                # create lookup table from tespy model
                self.lut_charge_path = self.wdir + self.sc + '_lut_charge.csv'
                self.lut_discharge_path = self.wdir + self.sc + '_lut_discharge.csv'
                self.tespy_create_lut()

            else:
                # use existing lut at specified path
                self.lut_charge_path = self.wdir + self.lut_charge_path
                self.lut_discharge_path = self.wdir + self.lut_discharge_path

        if self.method == 'lut':
            self.load_lut_model()

    def load_tespy_model(self):

        # load tespy models with the network_reader module
        self.tespy_charge = load_network(self.wdir + self.tespy_charge_path)
        self.tespy_charge.set_attr(iterinfo=False)
        self.tespy_charge.solve('design', init_only=True)
        self.tespy_discharge = load_network(self.wdir + self.tespy_discharge_path)
        self.tespy_discharge.set_attr(iterinfo=False)
        self.tespy_discharge.solve('design', init_only=True)

        self.power_plant_layout()

    def load_lut_model(self):
        # load splines from .csv data
        self.lut_charge = self.load_lookup_table(self.lut_charge_path)
        self.lut_discharge = self.load_lookup_table(self.lut_discharge_path)
        if (max(self.lut_charge[1]) != max(self.lut_discharge[1]) or
            min(self.lut_charge[1]) != min(self.lut_discharge[1])):
            msg = 'Pressure limits for charging and discharging do not match!'
            logging.error(msg)
            raise ValueError(msg)

        self.p_max = max(self.lut_charge[1])
        self.p_min = min(self.lut_charge[1])

    def power_plant_layout(self):
        """
        Power plant layout calculation to determine power plant design point using
        nominal power input/output and nominal pressure as inputs.
        """
        msg = 'Starting power plant layout calculation, compressor part.'
        logging.debug(msg)
        # charging
        self.tespy_charge.busses[self.power_bus_charge].set_attr(P=self.power_nominal_charge)
        self.tespy_charge.connections[self.pressure_conn_charge].set_attr(
                p=self.pressure_nominal_charge,
                m=ref(self.tespy_charge.connections[self.massflow_conn_charge], 1 / self.num_wells, 0))
        self.tespy_charge.connections[self.massflow_conn_charge].set_attr(m=np.nan)
        self.tespy_charge.components[self.pipe_charge].set_attr(L=self.min_well_depth)
        self.tespy_charge.solve('design')
        self.tespy_charge.save(self.wdir + self.sc + '_charge_design')
        self.m_nom_charge = self.tespy_charge.connections[self.massflow_conn_charge].m.val_SI
        self.m_min_charge = self.m_nom_charge * self.massflow_min_rel
        self.m_max_charge = self.m_nom_charge * self.massflow_max_rel
        msg = 'Nominal mass flow for charging is ' + str(self.m_nom_charge) + ' at nominal power ' + str(self.power_nominal_charge) + ' and nominal pressure ' + str(self.pressure_nominal_charge) + '.'
        logging.debug(msg)

        msg = 'Starting power plant layout calculation, turbine part.'
        logging.debug(msg)
        # discharging
        self.tespy_discharge.busses[self.power_bus_discharge].set_attr(P=self.power_nominal_discharge)
        self.tespy_discharge.connections[self.pressure_conn_discharge].set_attr(
                p=self.pressure_nominal_discharge,
                m=ref(self.tespy_discharge.connections[self.massflow_conn_discharge], 1 / self.num_wells, 0))
        self.tespy_discharge.connections[self.massflow_conn_discharge].set_attr(m=np.nan)
        self.tespy_charge.components[self.pipe_discharge].set_attr(L=self.min_well_depth)
        self.tespy_discharge.solve('design')
        self.tespy_discharge.save(self.wdir + self.sc + '_discharge_design')
        self.m_nom_discharge = self.tespy_discharge.connections[self.massflow_conn_discharge].m.val_SI
        self.m_min_discharge = self.m_nom_discharge * self.massflow_min_rel
        self.m_max_discharge = self.m_nom_discharge * self.massflow_max_rel
        msg = 'Nominal mass flow for discharging is ' + str(self.m_nom_discharge) + ' at nominal power ' + str(self.power_nominal_discharge) + ' and nominal pressure ' + str(self.pressure_nominal_discharge) + '.'
        logging.debug(msg)

    def tespy_create_lut(self):
        """
        Creates a lookup table for mass flow, pressure and power (charging and discharging).
        The lookup table is based on the power plant layout calculation.
        """
        # resolution of the lookup table data grid
        grid_num = 8

        # charging
        msg = 'Generating lookup table for charging from TESPy power plant model.'
        logging.debug(msg)

        m_range = np.linspace(self.m_min_charge, self.m_max_charge, grid_num)[::-1]
        pressure_range = np.linspace(self.p_min, self.p_max, grid_num)[::-1]

        df = pd.DataFrame(columns=pressure_range)
        df_heat = pd.DataFrame(columns=pressure_range)
        self.tespy_charge.busses[self.power_bus_charge].set_attr(P=np.nan)

        for m in m_range:
            power = []
            heat = []
            self.tespy_charge.connections[self.massflow_conn_charge].set_attr(m=m)
            for p in pressure_range:
                self.tespy_charge.connections[self.pressure_conn_charge].set_attr(p=p)
                self.tespy_charge.solve(mode='offdesign', design_path=self.wdir + self.sc + '_charge_design')
                if self.tespy_charge.res[-1] > 1e-3:
                    msg = 'Error on lookup table creation: Could not find a solution for input pair: mass flow=' + str(round(m, 1)) + ', pressure=' + str(round(p, 3)) + '.'
                    logging.error(msg)
                    raise TESPyNetworkError(msg)
                else:
                    power += [self.tespy_charge.busses[self.power_bus_charge].P.val]
                    heat += [self.tespy_charge.busses[self.heat_bus_charge].P.val]

            df.loc[m] = power
            df_heat.loc[m] = heat

        df.to_csv(self.wdir + self.sc + '_lut_charge_power.csv')
        df_heat.to_csv(self.wdir + self.sc + '_lut_charge_heat.csv')

        # discharging
        msg = 'Generating lookup table for charging from TESPy power plant model.'
        logging.debug(msg)

        m_range = np.linspace(self.m_min_discharge, self.m_max_discharge, grid_num)[::-1]
        pressure_range = np.linspace(self.p_min, self.p_max, grid_num)[::-1]

        df = pd.DataFrame(columns=pressure_range)
        df_heat = pd.DataFrame(columns=pressure_range)
        self.tespy_discharge.busses[self.power_bus_discharge].set_attr(P=np.nan)

        for m in m_range:
            power = []
            heat = []
            self.tespy_discharge.connections[self.massflow_conn_discharge].set_attr(m=m)
            for p in pressure_range:
                self.tespy_discharge.connections[self.pressure_conn_discharge].set_attr(p=p)
                self.tespy_discharge.solve(mode='offdesign', design_path=self.wdir + self.sc + '_discharge_design')
                if self.tespy_discharge.res[-1] > 1e-3:
                    msg = 'Error on lookup table creation: Could not find a solution for input pair: mass flow=' + str(round(m, 1)) + ', pressure=' + str(round(p, 3)) + '.'
                    logging.error(msg)
                    raise TESPyNetworkError(msg)
                else:
                    power += [self.tespy_discharge.busses[self.power_bus_discharge].P.val]
                    heat += [self.tespy_discharge.busses[self.heat_bus_discharge].P.val]

            df.loc[m] = power
            df_heat.loc[m] = heat

        df.to_csv(self.wdir + self.sc + '_lut_discharge_power.csv')
        df_heat.to_csv(self.wdir + self.sc + '_lut_discharge_heat.csv')

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
            msg = 'Pressure is below minimum pressure: min=' + str(self.p_min) + ', value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0, 0
        if pressure - 1e-4 > self.p_max:
            msg = 'Pressure is above maximum pressure: max=' + str(self.p_max) + ', value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0, 0

        if mode == 'shut-in':
            return 0, 0, 0

        if self.method == 'tespy':
            if mode == 'charging':
                # if power too small
                if abs(power) < abs(self.power_nominal_charge / 100):
                    return 0, 0, 0

                design_path = self.wdir + self.sc + '_charge_design'
                # set power of bus
                self.tespy_charge.busses[self.power_bus_charge].set_attr(P=power)
                # set pressure at interface
                self.tespy_charge.connections[self.pressure_conn_charge].set_attr(p=pressure)
                self.tespy_charge.connections[self.massflow_conn_charge].set_attr(m=np.nan)

                try:
                    self.tespy_charge.solve(mode='offdesign', design_path=design_path)
                    m = self.tespy_charge.connections[self.massflow_conn_charge].m.val_SI
                    heat = self.tespy_charge.busses[self.heat_bus_charge].P.val

                    if self.tespy_charge.res[-1] > 1e-3:
                        msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                        print(msg)
                        logging.error(msg)
                        return 0, 0, 0
                    elif m < self.m_min_charge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' below minimum mass flow.'
                        print(msg)
                        logging.error(msg)
                        return 0, 0, 0
                    elif m > self.m_max_charge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' above maximum mass flow. Adjusting power to match maximum allowed mass flow.'
                        print(msg)
                        logging.warning(msg)
                        return self.get_power(self.m_max_charge, pressure, mode)
                    else:
                        msg = 'Calculation successful for power=' + str(power) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                        print(msg)
                        logging.debug(msg)
                        return m, power, heat

                except:
                    # except general errors in calculation
                    msg = 'ERROR: Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                    print(msg)
                    logging.error(msg)
                    return 0, 0, 0

            elif mode == 'discharging':
                if abs(power) < abs(self.power_nominal_discharge / 100):
                    return 0, 0, 0

                design_path = self.wdir + self.sc + '_discharge_design'
                # set power of bus
                self.tespy_discharge.busses[self.power_bus_discharge].set_attr(P=power)
                # set pressure at interface
                self.tespy_discharge.connections[self.pressure_conn_discharge].set_attr(p=pressure)
                self.tespy_discharge.connections[self.massflow_conn_discharge].set_attr(m=np.nan)

                try:
                    self.tespy_discharge.solve(mode='offdesign', design_path=design_path, init_path=design_path)
                    m = self.tespy_discharge.connections[self.massflow_conn_discharge].m.val_SI
                    heat = self.tespy_discharge.busses[self.heat_bus_discharge].P.val

                    if self.tespy_discharge.res[-1] > 1e-3:
                        msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                        print(msg)
                        logging.error(msg)
                        return 0, 0, 0
                    elif m < self.m_min_discharge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' below minimum mass flow.'
                        print(msg)
                        logging.error(msg)
                        return 0, 0, 0
                    elif m > self.m_max_discharge:
                        msg = 'Mass flow for input pair power=' + str(power) + ' pressure=' + str(pressure) + ' above maximum mass flow. Adjusting power to match maximum allowed mass flow.'
                        print(msg)
                        logging.warning(msg)
                        return self.get_power(self.m_max_discharge, pressure, mode)
                    else:
                        msg = 'Calculation successful for power=' + str(power) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                        print(msg)
                        logging.debug(msg)
                        return m, power, heat

                except:
                    # except general errors in calculation
                    msg = 'Could not find a solution for input pair power=' + str(power) + ' pressure=' + str(pressure) + '.'
                    print(msg)
                    logging.error(msg)
                    return 0, 0, 0

            else:
                raise ValueError('Mode must be charging or discharging.')

        elif self.method == 'lut':
            if mode == 'charging':
                func = self.lut_charge
                f = 1

            elif mode == 'discharging':
                func = self.lut_discharge
                f = -1
                power = abs(power)
            else:
                raise ValueError('Mode must be charging or discharging.')

            # find position of pressure in pressure array (p_pos)
            # and calculate power vector pow_arr
            p_pos = np.searchsorted(func[1], pressure)
            if p_pos == len(func[1]):
                pow_arr = func[2][:, p_pos - 1]
            elif p_pos == 0:
                pow_arr = func[2][:, 0]
            else:
                pow_frac = (pressure - func[1][p_pos - 1]) / (func[1][p_pos] - func[1][p_pos - 1])
                pow_arr = func[2][:, p_pos - 1] + pow_frac * (func[2][:, p_pos] - func[2][:, p_pos - 1])

            # find position of power in power array (pow_pos)
            # and calculate mass flow
            pow_pos = np.searchsorted(pow_arr, power)
            if pow_pos == len(pow_arr):
                # mass flow is larger than last value in array
                # throttle power to maximum power available for given pressure value
                msg = 'Mass flow for input pair power=' + str(power * f) + ' pressure=' + str(pressure) + ' above maximum mass flow. Adjusting power to match maximum allowed mass flow.'
                logging.warning(msg)
                m = func[0][pow_pos - 1]
                msg = 'Calculation successful for power=' + str(max(pow_arr) * f) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                logging.debug(msg)
                return m, max(pow_arr) * f
            elif pow_pos == 0:
                # mass flow is identical or lower than minimum mass flow
                if power < min(pow_arr) * 0.999:
                    # shut plant down
                    msg = 'Could not find a solution for input pair power=' + str(power * f) + ' pressure=' + str(pressure) + '.'
                    logging.error(msg)
                    return 0, 0
                else:
                    # minimum mass flow
                    m = func[0][0]
                    msg = 'Calculation successful for power=' + str(power * f) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                    logging.debug(msg)
                    return m, power * f
            else:
                m_frac = (power - pow_arr[pow_pos - 1]) / (pow_arr[pow_pos] - pow_arr[pow_pos - 1])
                m = func[0][pow_pos - 1] + m_frac * (func[0][pow_pos] - func[0][pow_pos - 1])
                msg = 'Calculation successful for power=' + str(power * f) + ' pressure=' + str(pressure) + '. Mass flow=' + str(m) + '.'
                logging.debug(msg)
                return m, power * f
        else:
            raise ValueError('Method must be tespy or lut.')

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
        mass_flow_actual : float
            Actual mass flow of the power plant.

        power : float
            Actual electrical power input/output of the power plant.
        """
        if pressure + 1e-4 < self.p_min:
            msg = 'Pressure is below minimum pressure: min=' + str(self.p_min) + ', value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0, 0
        if pressure - 1e-4 > self.p_max:
            msg = 'Pressure is above maximum pressure: max=' + str(self.p_max) + ', value=' + str(pressure) + '.'
            logging.error(msg)
            return 0, 0, 0

        if mode == 'shut-in':
            return 0, 0, 0

        if self.method == 'tespy':
            if mode == 'charging':
                if mass_flow < self.m_min_charge - 1e-4:
                    msg = 'Mass flow is below minimum mass flow, shutting down power plant.'
                    print(msg)
                    logging.error(msg)
                    return 0, 0, 0
                elif mass_flow > self.m_max_charge + 1e-4:
                    msg = 'Mass flow above maximum mass flow. Adjusting mass flow to maximum allowed mass flow.'
                    print(msg)
                    logging.warning(msg)
                    return self.get_power(self.m_max_charge, pressure, mode)


                design_path = self.wdir + self.sc + '_charge_design'
                # set power of bus
                self.tespy_charge.busses[self.power_bus_charge].set_attr(P=np.nan)
                # set pressure at interface
                self.tespy_charge.connections[self.pressure_conn_charge].set_attr(p=pressure)
                self.tespy_charge.connections[self.massflow_conn_charge].set_attr(m=mass_flow)

                self.tespy_charge.solve(mode='offdesign', design_path=design_path)

                power = self.tespy_charge.busses[self.power_bus_charge].P.val
                heat = self.tespy_charge.busses[self.heat_bus_charge].P.val
                msg = 'Calculation successful for mass flow=' + str(mass_flow) + ' pressure=' + str(pressure) + '. Power=' + str(power) + '.'
                print(msg)
                logging.debug(msg)
                return mass_flow, power, heat

            elif mode == 'discharging':
                if mass_flow < self.m_min_discharge - 1e-4:
                    msg = 'Mass flow is below minimum mass flow, shutting down power plant.'
                    print(msg)
                    logging.error(msg)
                    return 0, 0, 0
                elif mass_flow > self.m_max_discharge + 1e-4:
                    msg = 'Mass flow above maximum mass flow. Adjusting mass flow to maximum allowed mass flow.'
                    print(msg)
                    logging.warning(msg)
                    return self.get_power(self.m_max_discharge, pressure, mode)


                design_path = self.wdir + self.sc + '_discharge_design'
                # set power of bus
                self.tespy_discharge.busses[self.power_bus_discharge].set_attr(P=np.nan)
                # set pressure at interface
                self.tespy_discharge.connections[self.pressure_conn_discharge].set_attr(p=pressure)
                self.tespy_discharge.connections[self.massflow_conn_discharge].set_attr(m=mass_flow)

                self.tespy_discharge.solve(mode='offdesign', design_path=design_path)

                power = self.tespy_discharge.busses[self.power_bus_discharge].P.val
                heat = self.tespy_discharge.busses[self.heat_bus_discharge].P.val
                msg = 'Calculation successful for mass flow=' + str(mass_flow) + ' pressure=' + str(pressure) + '. Power=' + str(power) + '.'
                print(msg)
                logging.debug(msg)

                return mass_flow, power, heat

            else:
                raise ValueError('Mode must be charge or discharge.')

        elif self.method == 'lut':
            if mode == 'charging':
                func = self.lut_charge
                f = 1

            elif mode == 'discharging':
                func = self.lut_discharge
                f = -1
            else:
                raise ValueError('Mode must be charging or discharging.')

            if mass_flow < min(func[0]) - 1e-4:
                msg = 'Mass flow is below minimum mass flow, shutting down power plant.'
                logging.error(msg)
                return 0, 0
            elif mass_flow > max(func[0]) + 1e-4:
                msg = 'Mass flow above maximum mass flow. Adjusting mass flow to maximum allowed mass flow.'
                logging.warning(msg)
                return self.get_power(max(func[0]), pressure, mode)

            # find position of pressure in pressure array (p_pos)
            # and calculate power vector pow_arr
            p_pos = np.searchsorted(func[1], pressure)
            if p_pos == len(func[1]):
                pow_arr = func[2][:, p_pos - 1]
            elif p_pos == 0:
                pow_arr = func[2][:, 0]
            else:
                pow_frac = (pressure - func[1][p_pos - 1]) / (func[1][p_pos] - func[1][p_pos - 1])
                pow_arr = func[2][:, p_pos - 1] + pow_frac * (func[2][:, p_pos] - func[2][:, p_pos - 1])

            m_pos = np.searchsorted(func[0], mass_flow)
            if m_pos == 0:
                # minimum mass flow (above maximum mass flow not possible)
                msg = 'Calculation successful for mass flow=' + str(mass_flow) + ' pressure=' + str(pressure) + '. Power=' + str(pow_arr[0] * f) + '.'
                logging.debug(msg)
                return mass_flow, pow_arr[0] * f
            else:
                m_frac = (mass_flow - func[0][m_pos - 1]) / (func[0][m_pos] - func[0][m_pos - 1])
                power = (pow_arr[m_pos - 1] + m_frac * (pow_arr[m_pos] - pow_arr[m_pos - 1])) * f
                msg = 'Calculation successful for mass flow=' + str(mass_flow) + ' pressure=' + str(pressure) + '. Power=' + str(power) + '.'
                logging.debug(msg)
                return mass_flow, power

        else:
            raise ValueError('Method must be tespy or lut.')

    def load_lookup_table(self, path):
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
        y[y < 0] *= -1

        x1 = df.index.get_values()  # mass flow
        if x1[0] > x1[-1]:
            x1 = x1[::-1]
            y = y[::-1]
        x2 = np.array(list(map(float, list(df))))  # pressure
        if x2[0] > x2[-1]:
            x2 = x2[::-1]
            y = y[:, ::-1]

        return (x1, x2, y)
