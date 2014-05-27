
	             ,/k.
	            /  ih,		*MonteGrappa*
	       ,-' ,  `:7b 			v1.0
	     _.-/   '  /b.`.4p,
	  --   ,    ,-' ^6x, `.'^=._

		   
     WELCOME TO THE FIRST VERSION OF MONTEGRAPPA !    

=======================================================

This file contains all the necessary information to compile and install
Montegrappa v1.0.    



The DEFAULT version of the code needs gsl libraries to be installed in /opt/local/lib and their headers installed in /opt/local/include.
If they are in different locations, you have to modify the second and
third line of Makefile accordingly
->METTIAMO UNA VARIABILE *****

To compile the DEFAULT version of the code, just type

      $ make

An executable called "montegrappa" will be created in the root directory
of the sources.

To build the MPI version,  you must type 

      $ make version=MPI 

If you have a pre-builded copy of MonteGrappa on your machine, the command

      $ montegrappa

will print a welcome message with the current version.


 

