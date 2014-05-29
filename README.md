
	             ,/k.
	            /  ih,		*MonteGrappa*
	       ,-' ,  `:7b 			v1.0
	     _.-/   '  /b.`.4p,
	  --   ,    ,-' ^6x, `.'^=._

		   
     WELCOME TO THE FIRST VERSION OF MONTEGRAPPA !    

=======================================================

This file contains all the necessary information to compile and install
Montegrappa v1.0.    

The package contains two programs:
montegrappa: the Monte Carlo engine
grappino: a tool to create the input files for montegrappa.


The DEFAULT version of the code needs gsl libraries to be installed in /opt/local/lib and their headers installed in /opt/local/include.
If they are in different locations, you have to modify the second and
third line of Makefile accordingly
->METTIAMO UNA VARIABILE *****

To compile the DEFAULT version of montegrappa, just type

      $ make

An executable called "montegrappa" will be created in the root directory
of the sources.

To build the MPI version,  you must type 

      $ make version=MPI 

Finally, you can compile the tool grappino 
      
      $ make grappino


If you have a pre-compiled copy of the montecarlo engine on your machine, the command

      $ montegrappa

without arguments will print a welcome/usage message with the current version.





 

