#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 15:17:46 2018

@author: witte
"""

import pandas as pd
import numpy as np
from scipy import interpolate


def __main__():
    """
    This is the main function containing the loop
    """

    power_plant = poly()

    p = 1
    P = 10

    mass_flow_plant(P, p, power_plant)


def mass_flow_plant(power, pressure, method):
    """
    returns the mass flow for given storage pressure and plant's power
    """
    if power == 0:
        return 0
    elif power > 0:
        return method.mass_flow(power, pressure, 'charge')
    else:
        return method.mass_flow(-power, pressure, 'discharge')


class poly:

    def __init__(self):

        self.coeff = {}
        self.coeff['charge'] = {}
        self.coeff['charge']['00'] = 44.97
        self.coeff['charge']['10'] = 2.394
        self.coeff['charge']['01'] = -1.568
        self.coeff['charge']['20'] = -0.0002792
        self.coeff['charge']['11'] = -0.0106
        self.coeff['charge']['02'] = 0.0128

        self.coeff['discharge'] = {}
        self.coeff['discharge']['00'] = 49.24
        self.coeff['discharge']['10'] = 3.265
        self.coeff['discharge']['01'] = -0.9576
        self.coeff['discharge']['20'] = 0.001451
        self.coeff['discharge']['11'] = -0.01171
        self.coeff['discharge']['02'] = 0.01161

    def mass_flow(self, power, pressure, state):

        if state == 'charge':
            coeff = self.coeff['charge']
        else:
            coeff = self.coeff['discharge']

        return (coeff['00'] + coeff['10'] * power + coeff['01'] * pressure +
                coeff['20'] * power ** 2 + coeff['11'] * pressure * power +
                coeff['02'] * pressure ** 2)


class spline:

    def __init__(self, path):

        self.charge = self.load_lookup(path)
        self.discharge = self.load_lookup(path)

    def mass_flow(self, power, pressure, state):

        if state == 'charge':
            return self.charge(power, pressure)
        else:
            return self.discharge(power, pressure)

    def load_lookup(self, path):

        df = pd.read_csv(path, index_col=0)

        x1 = df.index.get_values()  # pressure
        x2 = np.array(list(map(float, list(df))))  # power
        y = df.as_matrix()  # mass flow

        func = interpolate.RectBivariateSpline(x1, x2, y)
        return func


class tespy:

    def __init__(self, tespy_charge, tespy_discharge):

        


def reverse_2d(params, y):
    r"""
    reverse function for lookup table

    :param params: variable function parameters
    :type params: list
    :param y: functional value, so that :math:`x_2 -
              f\left(x_1, y \right) = 0`
    :type y: float
    :returns: residual value of the function :math:`x_2 -
              f\left(x_1, y \right)`
    """
    func, x1, x2 = params[0], params[1], params[2]
    return x2 - func.ev(x1, y)


def reverse_2d_deriv(params, y):
    r"""
    derivative of the reverse function for a lookup table

    :param params: variable function parameters
    :type params: list
    :param y: functional value, so that :math:`x_2 -
              f\left(x_1, y \right) = 0`
    :type y: float
    :returns: partial derivative :math:`\frac{\partial f}{\partial y}`
    """
    func, x1 = params[0], params[1]
    return - func.ev(x1, y, dy=1)


def newton(func, deriv, params, k, **kwargs):
    r"""
    find zero crossings of function func with 1-D newton algorithm,
    required for reverse functions of fluid mixtures

    :param func: function to find zero crossing in
    :type func: function
    :param deriv: derivative of the function
    :type deriv: function
    :param params: vector containing parameters for func
    :type params: list
    :param k: target value for function func
    :type k: numeric
    :returns: val (float) - val, so that func(params, val) = k

    **allowed keywords** in kwargs:

    - val0 (*numeric*) - starting value
    - valmin (*numeric*) - minimum value
    - valmax (*numeric*) - maximum value
    - imax (*numeric*) - maximum number of iterations

    .. math::

        x_{i+1} = x_{i} - \frac{f(x_{i})}{\frac{df}{dx}(x_{i})}\\
        f(x_{n}) \leq \epsilon, \; n < 10\\
        n: \text{number of iterations}
    """

    # default valaues
    val = kwargs.get('val0', 300)
    valmin = kwargs.get('valmin', 70)
    valmax = kwargs.get('valmax', 3000)
    imax = kwargs.get('imax', 10)

    # start newton loop
    res = 1
    i = 0
    while abs(res) >= err:
        # calculate function residual
        res = k - func(params, val)
        # calculate new value
        val += res / deriv(params, val)

        # check for value ranges
        if val < valmin:
            val = valmin
        if val > valmax:
            val = valmax
        i += 1

        if i > imax:
            print('Newton algorithm was not able to find a feasible'
                  'value for function ' + str(func) + '.')

            break

    return val
