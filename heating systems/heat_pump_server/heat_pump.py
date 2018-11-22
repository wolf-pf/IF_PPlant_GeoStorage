#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:32:17 2018

@author: witte
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 08:44:36 2018

@author: witte
"""

from tespy import cmp, con, nwk

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# %% network

nw = nwk.network(fluids=['water', 'NH3', 'air'],
                 T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
                 p_range=[0.1, 100], T_range=[1, 500], h_range=[15, 3500])

# %% components

# sources & sinks

c_in = cmp.source('coolant in')
cb = cmp.source('consumer back flow')
cf = cmp.sink('consumer feed flow')
amb = cmp.source('ambient air')
amb_out1 = cmp.sink('sink ambient 1')
amb_out2 = cmp.sink('sink ambient 2')
c_out = cmp.sink('coolant out')

# ambient air system
fan = cmp.compressor('pump')

# consumer system

cd = cmp.condenser('condenser')
dhp = cmp.pump('district heating pump')
cons = cmp.heat_exchanger_simple('consumer')

# evaporator system

ves = cmp.vessel('vessel')
dr = cmp.drum('drum')
ev = cmp.heat_exchanger('evaporator')
su = cmp.heat_exchanger('superheater')
erp = cmp.pump('evaporator reciculation pump')

# compressor-system

cp1 = cmp.compressor('compressor 1')

# %% connections

# consumer system

c_in_cd = con.connection(c_in, 'out1', cd, 'in1')

cb_dhp = con.connection(cb, 'out1', dhp, 'in1')
dhp_cd = con.connection(dhp, 'out1', cd, 'in2')
cd_cons = con.connection(cd, 'out2', cons, 'in1')
cons_cf = con.connection(cons, 'out1', cf, 'in1')

nw.add_conns(c_in_cd, cb_dhp, dhp_cd, cd_cons, cons_cf)

# connection condenser - evaporator system

cd_ves = con.connection(cd, 'out1', ves, 'in1')

nw.add_conns(cd_ves)

# evaporator system

ves_dr = con.connection(ves, 'out1', dr, 'in1')
dr_erp = con.connection(dr, 'out1', erp, 'in1')
erp_ev = con.connection(erp, 'out1', ev, 'in2')
ev_dr = con.connection(ev, 'out2', dr, 'in2')
dr_su = con.connection(dr, 'out2', su, 'in2')

nw.add_conns(ves_dr, dr_erp, erp_ev, ev_dr, dr_su)

amb_fan = con.connection(amb, 'out1', fan, 'in1')
fan_su = con.connection(fan, 'out1', su, 'in1')
su_ev = con.connection(su, 'out1', ev, 'in1')
ev_amb_out = con.connection(ev, 'out1', amb_out1, 'in1')

nw.add_conns(amb_fan, fan_su, su_ev, ev_amb_out)

# connection evaporator system - compressor system

su_cp1 = con.connection(su, 'out2', cp1, 'in1')

nw.add_conns(su_cp1)

# compressor-system

cp1_c_out = con.connection(cp1, 'out1', c_out, 'in1')

nw.add_conns(cp1_c_out)

# %% component parametrization

# condenser system

cd.set_attr(pr1=0.99, pr2=0.99, ttd_u=5)
dhp.set_attr(eta_s=0.8)
cons.set_attr(pr=0.99, offdesign=['zeta'])

# ambient air

fan.set_attr(eta_s=0.7, pr=1.005, mode='man')

# evaporator system

ves.set_attr(mode='man')
ev.set_attr(pr1=0.99, pr2=0.99, ttd_l=10,
            kA_char1='EVA_HOT', kA_char2='EVA_COLD',
            design=['pr1', 'ttd_l'], offdesign=['zeta1', 'kA'])
su.set_attr(pr1=0.99, pr2=0.99, ttd_u=5)
erp.set_attr(eta_s=0.8)

# compressor system

cp1.set_attr(eta_s=0.8, mode='man')

# %% connection parametrization

# condenser system

c_in_cd.set_attr(fluid={'air': 0, 'NH3': 1, 'water': 0})
cb_dhp.set_attr(T=25, p=10, fluid={'air': 0, 'NH3': 0, 'water': 1})
cd_cons.set_attr(T=50)
cons_cf.set_attr(h=con.ref(cb_dhp, 1, 0), p=con.ref(cb_dhp, 1, 0))

# evaporator system cold side

erp_ev.set_attr(m=con.ref(ves_dr, 4, 0), p0=5)
su_cp1.set_attr(p0=5, h0=1700)

# evaporator system hot side

amb_fan.set_attr(T=19, p=1, fluid={'air': 0, 'NH3': 0, 'water': 1})
fan_su.set_attr()
ev_amb_out.set_attr(T=12)

# compressor-system

cp1_c_out.set_attr(p=con.ref(c_in_cd, 1, 0), h=con.ref(c_in_cd, 1, 0))

# %% key paramter

#cons.set_attr(Q=-400e3)

heat_bus = con.bus('heat input', P=270e3)
heat_bus.add_comps({'c': su, 'char': -1}, {'c': ev, 'char': -1})
nw.add_busses(heat_bus)

# %% Calculation

nw.solve('design')
nw.save('heat_pump')
nw.print_results()

T_range = [19]
x = []
y = []
heat_in = pd.read_csv('LastRZ.csv')
Q_range = np.array([135e3])[::-1]

for T in T_range:
    amb_fan.set_attr(T=T)
    eps = []
    i = 0

    for Q in heat_in['0']:
        if Q >= 100:
            heat_bus.set_attr(P=Q*1e3)
            nw.set_printoptions(print_level='err')
            nw.solve('offdesign', init_file='heat_pump/results.csv',
                     design_file='heat_pump/results.csv')

            if len(nw.errors) > 1 or nw.lin_dep:

                print('error')

            else:

                x += [Q]
                y += [-cd.Q.val/ (cp1.P.val + fan.P.val + erp.P.val) - 1]

            i +=1
            print(i)
            if i > 50:
                break

plt.scatter(x, y)
plt.show()
