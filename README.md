	 ______ __  __ _____ _      _____ 
	|  ____|  \/  |_   _| |    |_   _|
	| |__  | \  / | | | | |      | |  
	|  __| | |\/| | | | | |      | |  
	| |____| |  | |_| |_| |____ _| |_ 
	|______|_|  |_|_____|______|_____|

emili
============

```emili``` is an algorithmic framework to instantiate several different 
optimization algorithms to solve hard combinatorial optimization problems. 
The EMILI framework has been designed to be used in automatic algorithm design
systems that use a grammar based method to generate algorithms. To this end, 
the EMILI framework uses a parser to instantiate algorithms at run time.
The framework is focused on the implementation of single solution SLS algorithms
and in particular hybrid stochastic local search algorithms.

This software is [open source](http://opensource.org/) and is distributed
under the terms of the
[BSD 2-Clause License](http://opensource.org/licenses/BSD-2-Clause).

If you use ```emili``` in your scientific work, please cite the works
below; if you do a derivative work, you should say it is derivative and also
cite the papers below:

 *  Federico Pagnozzi and Thomas Stützle. 
    **Automatic design of hybrid stochastic local search algorithms for 
    permutation flowshop problems.** 
    European journal of operational research, 276(2), 409-421, 2019.
    DOI: [10.1016/j.ejor.2019.01.018](https://doi.org/10.1016/j.ejor.2019.01.018)

#### License ####

Copyright (c) 2015, Federico Pagnozzi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

### Usage ###

EMILI INSTANCE_FILE_PATH PROBLEM < ALGORITHM DESCRIPTION > [-it|-ro time] [rnds seed] [ps]

### How to compile ###
Use cmake for the configuration and then make.

An example: 

*from the source dir

$ mkdir build

$ cd build

$ cmake ../

$ make

### How to start contributing ###

1) Create a subdir for your problem.
2) Create classes for your problem and solution ( extend the classes Problem and Solution).
3) Create classes for the algorithmic components and/or 
   algorithm templates you want to use to generate algorithms. 
4) Create a Builder to tell the framework how to load a problem instance 
   and how to instantiate the components you created.
5) Use the framework to instantiate algorithms to solve your problem!

In the third Chapter of my PhD thesis ([available here](https://dipot.ulb.ac.be/dspace/bitstream/2013/294557/4/thesis.pdf)) you can learn how the framework is structured
and how to add components and templates.
Check the examples in the template dir to have a quick idea of what to do.
Check the components implemented for the permutation flowshop problem 
to learn more (unfortunately, the code there may not be commented).
For any question, feel free to contact me to the address written below.

New source code created in subdirectories, like the template dir, 
needs to be added to the list of source files in the CMakeList.txt file.
See the commented code at line 5 of CMakeList.txt for an example.

Please, do not modify the main classes in any way.
The framework is made to be interoperable, flexible and modular.
Any change you may do to the main classe may break these properties.
If you discover a bug or want to propose an improvement, feel free
to contact me at federico.pagnozzi@ulb.ac.be
