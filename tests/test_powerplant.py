#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

@author: witte
"""

# %% imports

from coupled_simulation import cp
from coupled_simulation import pp

cd = cp.coupling_data('/home/witte/nextcloud/Documents/Hochschule/Dissertation/Kraftwerkssimulation/angus_model_coupling/testdata/testcase.main_ctrl.json')
cd.powerplant_path = cd.powerplant_path.replace('\\', '/')

test = pp.model(cd, 700, 12, 80, 40)

print(test.get_mass_flow(test.power_nominal_charge, test.pressure_nominal_charge, 'charging'))
print(test.get_mass_flow(test.power_nominal_discharge, test.pressure_nominal_discharge, 'discharging'))
print(test.get_mass_flow(test.power_nominal_charge, 40, 'charging'))
print(test.get_mass_flow(test.power_nominal_discharge, 80, 'discharging'))