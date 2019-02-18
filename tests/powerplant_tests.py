#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% imports

from tespy import nwk, cmp, con, hlp
from coupled_simulation import cp, pp
from nose.tools import eq_
import numpy as np
import shutil

# %% component tests


class powerplant_tests:

    def setup(self):
        self.cd = cp.coupling_data('/home/witte/nextcloud/Documents/Hochschule/Dissertation/Kraftwerkssimulation/angus_model_coupling/testdata/testcase.main_ctrl.json')
        self.cd.powerplant_path = cd.powerplant_path.replace('\\', '/')
        well_depth = 750
        num_wells = 10
        p_min = 40
        p_max = 80
        self.model = pp.model(self.cd, well_depth, num_wells, p_max, p_min)

    def test_mass_flow_charge(self):

        eq_(self.model.get_mass_flow(self.model.power_nominal_charge, self.model.pressure_nominal_charge, 'charging'), self.model.tespy_charge.imp_conns[self.model.tespy_charge.massflow_conn_charge].m.design)
#print(test.get_mass_flow(test.power_nominal_discharge, test.pressure_nominal_discharge, 'discharging'))
#print(test.get_mass_flow(test.power_nominal_charge, 40, 'charging'))
#print(test.get_mass_flow(test.power_nominal_discharge, 80, 'discharging'))