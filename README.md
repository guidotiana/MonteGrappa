===============================================

                     ,/k.
                    /  ih,              
               ,-' ,  `:7b                      
             _.-/   '  /b.`.4p,
          --   ,    ,-' ^6x, `.'^=._


                 WELCOME TO
               MONTEGRAPPA 1.3 

===============================================


This file will help you to compile and install Montegrappa v1.3 .
The package contains three programs:

	- "montegrappa", the main Monte Carlo engine 
	- "grappino", a tool to create the input files for montegrappa.
	- "mhistogram", a data analysis tool
	- "mgp2pdb", a simple converter from .pol files to .pdb

The code needs GSL libraries and the MPI environment to be correctly 
installed on your machine. This is a requirement to enable STEMPERING
and PARALLEL TEMPERING functions. 

Please check these dependencies before installing the software, or skip 
to the next section if you want a custom installation ( a basic version 
of Montegrappa can be compiled without GSL and MPI).

By default, GSL libraries should be installed in /opt/local/lib and their
headers installed in /opt/local/include.
If they are in different locations, you have to modify the compilation 
flags in the Makefile accordingly.

The whole distribution can be quickly compiled just typing in your terminal

	$ make all

The binary files "montegrappa", "montegrappa_mpi", "grappino", "mhistogram", 
and "mgp2pdb" will be created in the ./bin directory. 


OPENMP (experimental): To activate OpenMP, change the following lines 
in the Makefile

	#CC=gcc
	CC=gcc -fopenmp
	#CCMP=/usr/lib64/openmpi/bin/mpicc
	CCMP=/usr/lib64/openmpi/bin/mpicc -fopenmp

and set the variable

	$ export OMP_NUM_THREADS=N

where N is the number of threads you want to use. N=1 means
no OpenMP parallelization. Performances are system and parameters
dependent: always test your configuration before to run long simulations.  

================ CUSTOM INSTALLATION ====================================


(1) MonteGrappa Single Core, no STEMPERING (no PARALLEL TEMPERING)

	$ make cleanobj
	$ make

An executable called "montegrappa" will be created in the ./bin directory.

(2) MonteGrappa Single Core, with STEPERING (no PARALLEL TEMPERING)

    - SYSTEM REQUIREMENTS: GSL libraries

	$ make cleanobj	
	$ make version=STEMPERING

As well as the (1) version, an executable called "montegrappa" will be created in the ./bin directory.


(3) MonteGrappa Multi Core, with PARALLEL TEMPERING (no STEMPERING)

    - SYSTEM REQUIREMENTS: MPI environment


	$ make cleanobj
	$ make version=MPI

An executable called "montegrappa_mpi" will be created in the ./bin directory.


(3) grappino and mhistogram --------------------------------------------------
 
The tools can be compiled alone, with the commands

      $ make grappino
      $ make mhistogram
	  $ make mgp2pdb

The executables "grappino" and "mhistogram" will be created in the ./bin directory.
