# Partice In Cell method for low density plasma simulation in Hall-Thurster engine

Software able to reproduce using object-oreinted programming tool, dynamics of non-collisional, low density plasma. Ions dynamics is studied through Largangian equation
of momentum and interact with agglomerates of electrons that generate an electrostatic field, whose potential is used eventually to compute the elctric field thath moves ions.

**Installation** 

To compile the code using Make it's enough to type:
	$ make
while being in the Makefile folder.

	$ make build
and
	$ make all
produce the same result.

**Execution** 

To run the executable it's possible to use:
	$ make run
that is equivalent to type:
	$ ./pic2d
pic2d being the name of the executable.

To clean there are two commands:
	$ make clean
removes all the object files created during the compilation,
but keeps the executable and the output files;
	$ make distclean
calls clean and then removes also:
	- the executable
	- txt and csv output files that may be generated during execution
	- png and fig files that can be generated using the graph.m script

**NOTE THAT**

The software Matlab2018b (or equivalent) is necessary to run the
graph.m script and produce the plots, but to obtain all the
numerical results only a c++ compiler that supports the c++11 
standard it's strictly necessary. 

The code has been tested only with GNU Make and gcc, in particular using:
 -the APC course VM 
 -GNU Make 4.1, gcc 6.3 on Debian 9.


