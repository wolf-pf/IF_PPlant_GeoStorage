#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% imports

from tespy import nwk, cmp, con, hlp
from coupled_simulation import cp, pp
from nose.tools import eq_, raises
import numpy as np
import logging
import shutil

# %% component tests


class powerplant_tests:

    def setup(self):
        self.cd = cp.coupling_data('./testcase.main_ctrl.json')
        self.cd.powerplant_path = self.cd.powerplant_path.replace('\\', '/')
        well_depth = 750
        num_wells = 10
        p_min = 40
        p_max = 80
        self.model = pp.model(self.cd, well_depth, num_wells, p_max, p_min)

    def test_get_mass_flow_and_get_power(self):

        # mass flow in offdesign same as design mass flow?
        m_calc = round(self.model.get_mass_flow(self.model.power_nominal_charge, self.model.pressure_nominal_charge, 'charging')[0], 5)
        m_precalc = round(self.model.tespy_charge.imp_conns[self.model.massflow_conn_charge].m.design, 5)
        eq_(m_calc, m_precalc)

        # mass flow in offdesign same as design mass flow?
        m_calc = round(self.model.get_mass_flow(self.model.power_nominal_discharge, self.model.pressure_nominal_discharge, 'discharging')[0], 5)
        m_precalc = round(self.model.tespy_discharge.imp_conns[self.model.massflow_conn_discharge].m.design, 5)
        eq_(m_calc, m_precalc)

        # other

        shutil.rmtree('./charge_design', ignore_errors=True)
        shutil.rmtree('./discharge_design', ignore_errors=True)
