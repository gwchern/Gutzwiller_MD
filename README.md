# Gutzwiller_MD

This program is to perform Guzwiller molecular dynamics simulations on a toy model.

## compile
Download the GMD folder, go to the GMD folder, and type "make" to compile, the "GMD/Makefile" will take care of the rest. This program requires "armadillo" package installed, see http://arma.sourceforge.net/.

## run
Once the compilation is done, it will generate an executable file "gmd", type "./gmd" to run.

"./gmd" can be attached with some input parameters, e.g. "./gmd 5 1.5 ...",

1st parameter is the Hubbard U, 

2nd is the temperature 

3-5th are the Gutzwiller solver parameters, usually do not have to change them

6th the number of atoms

7th Wigner-Seitz radius, which determines the density of the system

8th total steps in the simulation

9th waiting steps for the system to equilibrate, measurements of the system will performed after it

10th restart_tag, determines the particles' starting positions and velocities

11-12th, 11th quench_tag, whether to quench the temperature, if it does, 12th is the quench temperature

These parameters can be found in the beginning of the main() function.

## contents
gqmd.cpp contains the main() function

analysis.hpp and analysis.cpp mainly defines the measurement functons on the molecular dynamics (MD) simulation.

tbmd.hpp and tbmd.cpp is mainly to implement the MD process.

potential.hpp, potential.cpp, gutzwiller.hpp and gutzwiller.cpp are the implementation of the Gutzwiller solver for the correlatated electron system.

unit.hpp defined the units in the simulation.










