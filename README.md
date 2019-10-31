	 ______ __  __ _____ _      _____ 
	|  ____|  \/  |_   _| |    |_   _|
	| |__  | \  / | | | | |      | |  
	|  __| | |\/| | | | | |      | |  
	| |____| |  | |_| |_| |____ _| |_ 
	|______|_|  |_|_____|______|_____|

commit : e5dffeb8ea6a4577658872eac6ae82376edd42e5
Usage:

EMILI INSTANCE_FILE_PATH PROBLEM <ALGORITHM DESCRIPTION> [-it|-ro time] [rnds seed] [ps]


Please, do not modify the main classes in any way.
The framework is made to be interoperable, flexible and modular.
Any change you may do to the main classe may break these properties.
If you discover a bug or want to propose an improvement, feel free
to contact me at federico.pagnozzi@ulb.ac.be
 
How to compile
Use cmake for the configuration and then make.

An example: 

*from the source dir
# mkdir build
# cd build
# cmake ../
# make

How to proceed: 

1) Create a subdir for your problem.
2) Create classes for your problem and solution ( extend the classes Problem and Solution).
3) Create classes for the algorithmic components and/or 
   algorithm templates you want to use to generate algorithms. 
4) Create a Builder to tell the framework how to load a problem instance 
   and how to instantiate the components you created.
5) Use the framework to instantiate algorithms to solve your problem!

Use the manual (Chapter 3 of my thesis) to learn how the framework is structured
and how to add components and templates.
Check the examples in the template dir to have a quick idea of what to do.
Check the components implemented for the permutation flowshop problem 
to learn more (unfortunately, the code there may not be commented).
For any question, feel free to contact me to the address written above.

New source code created in subdirectories, like the template dir, 
needs to be added to the list of source files in the CMakeList.txt file.
See the commented code at line 5 of CMakeList.txt for an example.
