=======================================
0. Summary
=======================================

This program is a parallel island model genetic algorithm implementation, created for an undergrad senior project researching the parallelization and optimization of genetic algorithms for high-performance distributed computing systems. It was written in C and makes use of the OpenMPI library. That being said, this program is more functional than it is pretty.

It includes two simple fitness functions (for performance testing purposes only; these can obviously be solved by more direct means):  
Simple: Maximizes the function y^2 + x^2 - a^2 + b^2 - 2xy - 2ab  
Shortest 3-dimensional path: Finds the shortest path in 3-dimensional space between two specified points around a set of obstacles.

=======================================
1. Requirements
=======================================

This program requires OpenMPI: http://www.open-mpi.org

=======================================
2. Usage
=======================================

Compile and run ga from the command line. Simulation parameters can be entered either from user input or from previously saved parameter files (samples provided).

Configuration for reporting options and default GA parameters are located in config.h.

=======================================
3. User-specified values
=======================================

**Population size:** 
A positive integer value for the size of the population. For parallel execution, this determines the size of each sub-population.

**Fitness function:** 
The type of fitness function to use.   
Possible values:  
0: Maximizes a simple function  
1: Finds the shortest 3D path between 2 points that does not collide with a set of objects

**Termination type:** 
The condition on which the program will terminate.  
Possible values:  
0: Fixed number of generations  
1: Maximum fitness threshhold  
2: Average fitness threshhold  
3: Local convergence 

**End generation:** 
A positive integer value for the number of generations, if a termination type of fixed generations was selected.

**Fitness threshhold:** 
A positive value for an accepted fitness level, if a termination type of average or maximum fitness threshhold was selected.

**Convergence generation threshhold:** 
A positive integer value for the number of generations with low variation that should be considered a local convergence, if a termination type of local convergence was selected.

**Convergence variation:** 
A positive value for the variation allowed for a local convergence, if a termination type of local 
convergence was selected.

**RNG seed:** 
A value used to initialize the random rumber generator. 0 may be used to generate a new seed value. Mutliple runs of the program with the same seed value will yield the same results.

**Start and end points:** 
These define the end points between which a 3-dimensional path is contructed, entered in the format (x,y,z). Coordinates must be within SP_BOUND as defined in config.h.

**Number of obstacles:** 
These are the obstacles around which a path must be constructed. Each object is aproximated by a sphere, with a center entered in the format (x,y,z) and a positive integer value for the radius.


=======================================
4. Code Overview
=======================================

Code is laid out as follows:

**config.h:** 
Configurable parameters and program defaults  
**types.h:** 
Data structures used for this program  
**ga.c:** 
Main loop of program and genetic operators  
**fitness.c:** 
Fitness functions and associated helper methods  
**init.c:** 
Initialization and validation of starting populations  
**mt_mpi:** 
Parallel implementation of the Mersenne Twister RNG algorithm  
**report.c:** 
Helper functions for reporting population and fitness stats  

=======================================
5. References
=======================================

The following resources were used in the development of this program and its corresponding research paper:

Adeli, Hojjat. Cost Optimization of Structures: Fuzzy Logic, Genetic Algorithms, and Parallel Computing. Chichester: Wiley, 2006.
Cantu´-Paz, Erick. Efficient and Accurate Parallel Genetic Algorithms. Boston: Kluwer Academic, 2000.
Cantu´-Paz, Erick. Genetic and Evolutionary Computation: Genetic and Evolutionary Computation Conference: Proceedings. Berlin: Springer, 2003.
Cohoon, J.P., S.U. Hedge, W.N. Martin, D. Richards. Genetic Algorithms and Punctuated Equilibria in VLSI. Berlin: Springer-Verlag, 1991.
Danková, Martina. Approximation of Extensional Fuzzy Relations Over Residuated Lattices. University of Ostrava Institute for Research and Applications of Fuzzy Modeling: Ostrava, 2008.
De Jong, Kenneth Alan. An Analysis of the Behavior of a Class of Genetic Adaptive Systems. University of Michigan Press: Ann Arbor, 1975.
Falkenauer, Emanuel. Genetic Algorithms and Grouping Problems. Chichester: Wiley, 1998.
Gwiazda, Tomasz D. Genetic Algorithms Reference, Volume 1: Crossover for single-objective numerical optimization problems. Boston: Twayne Publishers, 2006. 
Haupt, Randy L. and Sue Ellen. Practical Genetic Algorithms, Second Edition. Hoboken: Wiley, 2004.
Schaeferm, Robert. Foundation of Global Genetic Optimization. Berlin: Springer, 2007.
Stender, Joachim. Parallel Genetic Algorithms: Theory & Applications. Amsterdam: IOS, 1993.
Vose, Michael D. The Simple Genetic Algorithm: Foundations and Theory. Cambridge: MIT Press, 1999.
