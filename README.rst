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
As TESPy is a equation-based approach, the system of equations is influenced by parameter specification of the user, too.



Plant layout
^^^^^^^^^^^^

Offdesign operation
^^^^^^^^^^^^^^^^^^^

Proxy model
+++++++++++

The main reason for implementing a proxy model in the interface is reduction of calculation time. The proxy model is a two dimensional lookup table,
linking the key figures of the compressed air energy storage to each other.

Software tests
--------------
