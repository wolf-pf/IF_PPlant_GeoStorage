~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Documentation of the Power Plant Geostorage Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
    :depth: 1
    :local:
    :backlinks: top

With the Power Plant Geostorage Interface it is possible to couple the simulation of storage and power plant operation.
The interface provides a data exchange between the sotrage and the power plant model and allows to utilize different simulations models.

Coupled Simulation
------------------

Description of the interface
++++++++++++++++++++++++++++

- modelinitialisation (read input time series, control files, prepare models)
- start loop for timestep

	- read power input/output from input time series, initial pressure from geological model
	- calculate input/output mass flow and actual power with power plant model (at given pressure level and power requirement)
	- calculate actual pressure and mass flow with geological storage model
	- check convergence:

		- assumed power of power plant model vs. actual pressure (go back to 2.2, recalculate mass flow)
		- pressure limits of storage (go back to 2.2, reduce mass flow, calculate power input/output)

	- write data to output file, if convergence aachieved or number of iterations exceeded:

		- timestamp
		- power
		- pressure
		- mass flow rate
		- temperature at boreholen (optional)


Input data
++++++++++

- general configuration for coupled simulation in scenarioname.main_ctrl.json
- time series containing planned power plant dispatch (input/output power)
- configuration of the power plant model (TESPy or proxy model) in powerplant/scenarioname.powerplant_ctrl.json

	- power plant setup
	- maximum/minimum mass flow
	- maximum/minimum pressure
	- path to power plant models
	- ...

- configuration of the geological storage model (ECLIPSE or proxy model) in geostorage/scenarioname.geostorage_ctrl.json

	- number of wells
	- geological setting
	- storage preparation
	- ...

Output data
+++++++++++

- output time series

	- timestamp
	- power
	- pressure
	- mass flow rate
	- temperature at boreholen (optional)

Geological Storage Module
-------------------------

Power Plant Module
------------------

The power plant module handles the power plant simulation. There are two main functions of the power plant module:

- Calculate the mass flow pressed into or extracted from the storage for a given pressure level at the bottom of the borehole and given power input or output respectively.
- Calculate the power input (or output) at given mass flow and given pressure level. This function is e. g. required, if the scheduled power can not be reached due to restrictions of the storage or the power plant.

For this purpose a power plant model is based on single components, e. g. compressors, turbines (turbine stages), heat exchangers, control valves or piping, is set up.
The interface implementation supports two different approaches for calculation of interface functionalities as mentioned above:

- Importing a TESPy power plant model or
- using a two dimensional lookup table linking the three quantities (mass flow, pressure level and power) as well as an inverse function.

The following sections will introduce the different modeling approaches and outline their strenghts and weaknesses.

TESPy model
+++++++++++

TESPy (Thermal Engineering Systems in Python) is a modular and component based power plant design and simulation software.
The Software is designed to simulate stationary operation of power plants, district heating systems, heat pumps or similar applications.

By connecting different components with each other a network is created, which can be represented by a system nonlinear of equations.
As TESPy is a equation-based approach, the system of equations is additionally influenced by parameter specification of the user.
In order to calculate the stationary operation of the plant, the fluid properties and the fluid composition of all connections between the components of the network has to be identified.
Thus, the system's variables are: mass flow, pressure, enthalpy and mass fraction of the fluid components. Component parameters can optionally be used as additional variables of the system, e. g. for the layout of a pipe diameter.
In order to numerically solve the system of equations, TESPy uses the multi dimensonal Newton–Raphson method.

The value of the system variables is calculated accoring to equation (0.0) in every iteration i:

.. math::

	\vec{x}_{i+1}=\vec{x}_i-J\left(\vec{x}_i\right)^{-1}\cdot f\left(\vec{x}_i\right)

Therefor, the calculation of the residual values of the equations :math:`f\left(\vec{x}_i\right)` as well as the calculation of the inversed jacobian matrix :math:`J\left(\vec{x}_i\right)` is required.
The algorithm is terminated, if the magnitude of the equations (vector norm :math:`||f\left(\vec{x}_i\right)||`) is smaller than a specified residual value:

.. math::

	||f(\vec{x}_i)|| > \epsilon
	
[1] provides further detailed information on TESPy.
	
For the interface, the power plant model is loaded with the tespy.network_reader-module, allowing to load the plant's topology and parametrisation as tabular data (from .csv-format).
After loading the plant model it is still possible to change the following parameters:

- depth of the wells :math:`L_{wells}` and number of wells :math:`n_{wells}`, as well as minimum and maximum pressure at the bottom of bore holes :math:`p_{min}`, :math:`p_{max}` provided by the geostorage model controle files.
- the nominal power, nominal pressure at bottom of the bore hole, maximum and minimum (relative) mass flow (in regard to mass flow at nominal power and pressure) provided by the powerplant model control files.

Based on these settings a plant design layout will be performed. All further operation will then reference this design point.
After the plant's design, the bottom bore hole pressure, the mass flow and the total power input/output are be the only exchangable parameters.
As mentioned in the introducing part, two different ways to control the power plant operation are required, which are outlined in the following sections.

Calculation of mass flow
^^^^^^^^^^^^^^^^^^^^^^^^

The determination of the mass flow rate represents the common operation mode, if the operation schedule of the power plant does not interfere with any restrictions of the geological storage.
The mass flow is calculated as a function of borehole pressure and electrical power: The electrical power input of the motor (power output of generator for discharging mode) is set to the target value from the input time series.
The pressure level at the bottom bore hole is retrieved by the geological storage simulation. Also, if the mass flow rate was specified in a prior calculation it will be unset for this case.
Following, TESPy will solve the plant model, whereby different outcomes are possible (see table ...).

The expected result is, that the TESPy solver is able to find a feasible solution for plant's point of operation. The calculated mass flow rate and the scheduled power are returned in this first case.
All other cases represent errors in the calculation or violations of the plant's operation limits. If the mass flow rate is higher than the allowed maximum mass flow rate (according to the power plant specifications), the power plant
will reduce its power to the value according to the maximum mass flow rate. Thus, maximum mass flow rate and a corrected value for the power are returned.

Other possible cases are pressure limit violations (higher than maximum or lower than minimum pressure), solver unable to find feasible solution or other errors within the calculation process.
In these cases, the power plant is shut down and a mass flow rate of 0 kg/s along with a power of 0 W is returned.

case;description;returned mass flow;returned power
calculation successful;TESPy solver found feasible solution;mass flow rate;scheduled power
mass flow > maximum mass flow;mass flow higher than possible;maximum mass flow;calculate power(maximum mass flow, pressure)
mass flow < minimum mass flow;mass flow too low for plant operation;0;0
p < p_min or p > p_max;pressure out of pressure limits;0;0
residual > 1e-3;no feasible solution for steady state was found by solver;0;0
other errors in calculation;Fluid property errors, unexpected errors in the solution process, other;0;0

Calculation of power
^^^^^^^^^^^^^^^^^^^^

The calculation of the power at a given mass flow rate and pressure is required, if the operation schedule of the storage leads to restrictions in mass flow rate.
If a mass flow rate calculated in a prior iteration can not be met due to pressure limitations of the geological storage or limitations of the power plant, the
actual power will be calculated in order show the deviation from the target.

For the calculation of the electrical power, mass flow rate and pressure are specified in the TESPy model. Possible errors in the calculation are identical to the
errors in the calculation of the mass flow rate (see table ...). In case of a successful calculation the calculated electrical power according to given mass flow rate
and pressure is returned.

Proxy model
+++++++++++

The main reason for implementing a proxy model in the interface is reduction of calculation time. The proxy model is a two dimensional lookup table,
linking the key figures of the compressed air energy storage to each other.

Software tests
--------------

Literature
----------
Francesco Witte. (2019, February 2). Thermal Engineering Systems in Python (Version latest). Zenodo.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2555867.svg
   :target: https://doi.org/10.5281/zenodo.2555867
