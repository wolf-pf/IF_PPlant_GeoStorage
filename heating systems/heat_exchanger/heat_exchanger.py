#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 08:44:36 2018

@author: witte
"""

from tespy import cmp, con, nwk

import numpy as np
import pandas as pd
from time import sleep
from matplotlib import pyplot as plt

heat_input = pd.read_csv('Zeitreihe Speicher.csv', decimal=',', sep=';')

# %% network

nw = nwk.network(fluids=['water'],
                 T_unit='C', p_unit='bar', h_unit='kJ / kg', m_unit='kg / s',
                 p_range=[0.1, 100], T_range=[1, 500], h_range=[15, 3500])

# %% components

# sources & sinks

tesin = cmp.sink('TES in')
tesout = cmp.source('TES out')
hsin = cmp.sink('HS in')
hsout = cmp.source('HS out')
he = cmp.heat_exchanger('heat exchanger')

# %% connections

# consumer system

tes_he = con.connection(tesout, 'out1', he, 'in2')
he_tes = con.connection(he, 'out2', tesin, 'in1')

hs_he = con.connection(hsout, 'out1', he, 'in1')
he_hs = con.connection(he, 'out1', hsin, 'in1')

nw.add_conns(tes_he, he_tes, hs_he, he_hs)

he.set_attr(pr1=0.98, pr2=0.98, ttd_u=5,
            design=['pr1', 'pr2', 'ttd_u'],
            offdesign=['zeta1', 'zeta2', 'kA'])

hs_he.set_attr(T=90, p=5, fluid={'water': 1})
he_hs.set_attr(T=60)

tes_he.set_attr(p=5, fluid={'water': 1})

# %% Schnittstellenparameter

tes_he.set_attr(T=40)

Q0 = 80e3
he.set_attr(Q=-Q0)

# %% Calculation

nw.solve('design')
nw.save('test')
nw.print_results()

he.set_attr(Q=-7e3)
nw.solve('offdesign', design_file='test/results.csv')
nw.save('low_Q')

m = []
q = []

for Q in heat_input['charge_Q']:

    if Q > 1e-1:
        he.set_attr(Q=-Q * 1e3)
        nw.set_printoptions(print_level='err')
        if he.Q.val > -10e3:
            nw.solve('offdesign', init_file='low_Q/results.csv', design_file='test/results.csv')
        else:
            nw.solve('offdesign', init_file='test/results.csv', design_file='test/results.csv')
#            if nw.lin_dep:
        print(round(Q, 3), round(tes_he.m.val, 3))

        m += [tes_he.m.val_SI]
        q += [Q]

plt.plot(q, m, 'x')
plt.show()


