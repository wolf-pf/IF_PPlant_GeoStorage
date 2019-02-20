~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Documentation of the Power Plant Geostorage Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. contents::
    :depth: 2
    :local:
    :backlinks: top
	
Abstract
--------
	
Introduction
------------

With the Power Plant Geostorage Interface it is possible to couple the simulation of storage and power plant operation.
The interface provides a control logic for the energy storage operation and exchanges data between the sotrage and the power plant model.
It also allows the exchange of the simulation models.

Coupling Module
---------------

Description of the interface
++++++++++++++++++++++++++++

The coupled simulation of the power plant and the geological storage is performed by exchanging the physical parameters mass flow and pressure between individual models.
The operation schedule is provided as input time series, thus, the models are coupled time stepwise and iterate in a loop until the interface parameters - bottom bore hole pressure and mass flow - are matched in both models:

The power plant model calculates the mass flow into or from the geological storage required in order to meet the scheduled power input or output given the actual bottom bore hole pressure.
Subsequently, the storage model calculates the bottom bore hole pressure for that mass flow and passes the pressure to the power plantmodel. This process is repeated until the models converge.

Convergence is defined by the following criterions:

.. math::

	\text{and}
	\begin{cases}
	\text{or} &
	\begin{cases}
	\frac{|p_{i-1} - p_{i}|}{p_{i}} < \epsilon \\
	|p_{i-1} - p_{i}| < \delta
	\end{cases} \\
	& \frac{|m_{pp} - m_{sto}|}{m_{sto}} < \epsilon
	\end{cases}	

If convergence is achieved, pressure, mass flow and actual power are written into the output file and the simulation for the next timestep is started.
		
Operation control
+++++++++++++++++

The operation schedule of the compressed air energy storage is an external input parameter for the model and not influenced by the coupled simulation.
But both, the geological storage and the power plant, have physical restrictions once they have been constructed.
These are e. g. pressure limits of the porous media or limitations in the power plant's range of operation (partload/overload).
If any of these limits are violated, the storage operation needs to be adjusted. This section outlines possible control

After processing the input parameters both models are initialised and the simulation starts. 

The interface 

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

The power plant module handles the power plant simulation in the interface. There are two main functions of the power plant module:

- Calculate the mass flow pressed into or extracted from the storage for a given pressure level at the bottom of the borehole and given power input or output respectively.
- Calculate the power input (or output) at given mass flow and given pressure level. This function is e. g. required, if the scheduled power can not be reached due to restrictions of the compressed air energy storage.

As mentioned in the first chapter, the basic task for the interface is to exchange data (pressure and mass flow) between a power plant and a storage model.
A very flexible representation of a power plant can be provided with two-dimensional lookup tables linking mass flow and pressure to power, as this approach
provides a standardised data structure for the power plant representation. The creation of tabular can be outsourced to any power plant simulation software.
Also, the calculation speed is much higher compared to running a power plant simulation software for each iteration. The downside of this approach is,
that the power plant design needs to be performed prior to the simulation and is not coupled directly to parameters of the geological storage. For example, these could be pressure limits, the number of wells or the depth of the wells.
Thus, the power plant module provides a second implementation of the power plant model for interface using the power plant simulation software TESPy (Thermal Engineering Systems in Python).
In this way, it is possible to design the power plant based on information of the geological storage. Additionally, the generation of tabular data is performed automatically based on these information.

The following two sections will describe the implementation of the tabular data power plant model and a component based power plant model using TESPy.

Tabular data model
++++++++++++++++++

<---
The main reason for implementing a proxy model in the interface is reduction of calculation time. The proxy model is a two dimensional lookup table,
linking the key figures of the compressed air energy storage to each other.
--->

Power plant simulation model
++++++++++++++++++++++++++++

TESPy is a modular and component based power plant design and simulation software.
The Software is designed to simulate stationary operation of power plants, district heating systems, heat pumps or similar applications.
By connecting different components with each other a network is created, which can be represented by a system nonlinear of equations.
	
In order to calculate the stationary operation of the system, the fluid properties and the fluid composition of all connections between the components of the network have to be identified.
Thus, the system's variables are: mass flow, pressure, enthalpy and mass fractions of the fluid components. Component parameters can optionally be used as additional variables of the system, e. g. for the layout of a pipe diameter.
In order to numerically solve the system of equations, TESPy uses the multi dimensonal Newton–Raphson method.

The value of the system variables is calculated accoring to equation (0.0) in every iteration i:

.. math::

	\vec{x}_{i+1}=\vec{x}_i-J\left(\vec{x}_i\right)^{-1}\cdot f\left(\vec{x}_i\right)

Therefore, the calculation of the residual values of the equations :math:`f\left(\vec{x}_i\right)` as well as the calculation of the inversed jacobian matrix :math:`J\left(\vec{x}_i\right)` is required.
The algorithm is terminated, if the magnitude of the equations (vector norm :math:`||f\left(\vec{x}_i\right)||`) is smaller than a specified residual value (eq. 0.0):

.. math::

	||f(\vec{x}_i)|| > \epsilon
	
Every component delivers a set of basic equations to the system of equations. Depending on the parametrisation of the components and connections more equations are added to the systems.
In this way, the topology and the parametrisation determine the set of equations used to describe the system: Different parametrisation of the same topological model results in a different system of equations.
Thus, we will not provide a detailed description of the power plant model here. Detailed information on the implementation of different components is provided in the online documentation and the API-documentation [1].

For usage in the interface the power plant model has to be designed by defining the topology and process parameters. The tespy.network_reader-module allows to load the plant's representation as tabular data (from .csv-format).
After loading the plant model it is still possible to change the following parameters:

- depth of the wells :math:`L_{wells}` and number of wells :math:`n_{wells}`, as well as minimum and maximum pressure at the bottom of bore holes :math:`p_{min}`, :math:`p_{max}` provided by the geostorage model controle files.
- the nominal power, nominal pressure at bottom of the bore hole, maximum and minimum (relative) mass flow (in regard to mass flow at nominal power and pressure) provided by the powerplant model control files.

The depth of the wells and the number of wells is used to determine the pressure losses in the pipes connecting the power plant's outlet (compression)/inlet (expansion) with the geological storage:
The total mass flow will be split up evenly amongst the pipes and the pressure loss is calculated using the darcy friction factor (equation 0.0).
Nominal power and nominal pressure are required to calculate the nominal mass flow.

.. math::

	0 = \delta p - \frac{\rho \cdot \lambda \left(Re, k_{s}, D \right) \cdot L \cdot c^2}{2 \cdot D}	

After the power plant design calculation, only the interface parameters (bottom bore hole pressure, total mass flow and total power input/output) are exchangable parameters.
The following two section describe the different ways to control the power plant operation using the interface parameters.

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

Software tests
--------------

Stuff
-----

<----
Each type of component represents a special case of an open thermodynamic system with a variable amount of inlets and outlets.
Equation 0.0 describes its energy balance. As the software solves for stationary operation, the sum of all mass flows from and into a component must be equal to zero (equation 0.0).

.. math::

	0 = \sum_i (\dot{m}_{out,i} \cdot h_{out,i}) - \sum_j (\dot{m}_{in,j} \cdot h_{in,j}) - \dot{Q} - P 

	0 = \sum_i \dot{m}_{out,i} - \sum_j \dot{m}_{in,j} \cdot h_{in,j}
--->

Literature
----------
Francesco Witte. (2019, February 2). Thermal Engineering Systems in Python (Version latest). Zenodo.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2555867.svg
   :target: https://doi.org/10.5281/zenodo.2555867
