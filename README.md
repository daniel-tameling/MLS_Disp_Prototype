MLS_Disp_Prototype
==================

prototype of the multilevel summation method for dispersion interactions

####################
## The compilation was tested with gnu and intel compiler
Everything necessary to run the prototype can be done via the Makefile
In order to compile the code there are two options: make  and  make noforce

make: after each time step, the output of the forces are written into a file, 
overwriting the previous time step.

I/O operation are especially costly, and increase the time spent in
the part of the code that also contains the integrator etc.  
To avoid I/O operations, compile with  make noforce

For the intel compiler some optimization flags can be added:
make intel  and  make noforce_intel

To run the program:
./msm_disp <file with particles> <config file> <lammps script used to create the particle-file>

The <file with particles> and its corresponding <lammps script used to create the particle-file> are located in ./input

The necessary MLS parameters have to be provided through the <config file>. 
Some example config files are in this folder.

HOWEVER, it is recommended to use the default options provided in the Makefile:
make run  and  make run500  and  make run32  and  make run256 and  make run_sh

make run500
1000 time steps 500-particle system

make run
1000 time steps 4,000-particle system

make run32
1000 time steps 32,000-particle system

make run256
1000 time steps for a 265,000-particle system

make run_sh 
compiles the code and runs a shell script that creates a plot.txt 
containing the potential energies (total, local, grid_based,
and self-interaction), and the force error for different cutoffs.