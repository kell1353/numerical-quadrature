# Project Overview
The program takes an individual real and integer iputs (chosen by the user) and calculates the neutron flux using the booles and monte carlo integration 
techniques. It calculates the flux in two different scerarios:

The first one, the flux is calculated over a solid box reactor where the distance of the detector is increasing away in distance from the reactor. The flux is then calculated
at each distance from reactor. The program also calculates the large `x_0` approximation of the neutron flux.

The second one, the flux is calculated over a hollow box reactor where the hollow space inside the reactor is the radius of the hollow sphere. The flux is then calculated
at each individual radius of the sphere between the user provided range for the radius. 

It calculates the results of the Boooles and Monte Carlo approximations for both scenarios and writes each scenarios results into seperate .dat to be evaluated in Jupyter.


# Compilation Instructions

Compilation is all done using the makefile in the repository. Type `make` into your command line to compile the files.

Initially sets nuclear-reactor as the executable name.
- Then creates the types object file.
- Then creates the quadrature object file using the types object file.
- Then creates the neutron_flux object file using the types, and quadrature object files.
- Then creates the read_write object file using the types, quadrature, and neutron_flux object files.
- Then creates the main object file using the types, and read_write object files.


# Usage Instructions

Once you have compiled everything execute the program (./nuclear-reactor)
The program will prompt you to type a real and integer inputs for various parameters.
It then will write the results of the program into a results_basic.dat file and a results_advanced.dat file. 
Once created you will be able to run the Jupyter file on this dataset. 

# Expected Behavior

Once you have typed your specified values in, it then will perform calculation and write the results of the program into a results_basic.dat file.
It will then prompt you to insert more values for the hollow box calculation and then it will write the results into the results_advanced.dat file.

The results_basic file should have 5 columns. `x_0`, 'booles', 'large `x_0`', 'monte carlo', 'MC uncertainty'.
- `x_0`: This column will contain the x coordinate of the detector's position.
- booles: This column will contain the booles quadrature integral of the neutron flux for each x_step.
- large `x_0`: This column will contain the large `x_0` approximation of the neutron flux.
- monte carlo: This column will contain the monte carlo integral evaluation of the neutron flux for each x_step.
- MC uncertainty: Provides the estimate of the uncertainty in the Monte Carlo integral.

The results_advanced file should have 5 columns. 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'.
- radius: This column will contain the radius of the hollow sphere.
- box booles: This column will contain the booles quadrature integral of the neutron flux for each x_step.
- hollow booles: This column will contain the difference of the box booles and sphere booles quadrature integrals for neutron flux for each radius.
- hollow monte carlo: This column will contain the hollow box monte carlo quadrature evaluation of the neutron flux.
- MC uncertainty: Provides the estimate of the uncertainty in the Monte Carlo integral.
