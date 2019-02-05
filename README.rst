~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Documentation of the Power Plant Geostorage Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

TESPy model
^^^^^^^^^^^

Proxy model
^^^^^^^^^^^

Software tests
--------------
