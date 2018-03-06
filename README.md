Biased Random-Key Genetic Algorithms for the Winner Determination Problem in 
Combinatorial Auctions
===============================

This project implements a biased random-key genetic algorithm (BRKGA) 
for the Winner Determination Problem in Combinatorial Auctions. The
algorithm details can be found in 

> C.E. Andrade, F.K. Miyazawa, M.C.G Resende, R.F. Toso. Biased Random-Key
> Genetic Algorithms for the Winner Determination Problem in Combinatorial
> Auctions. Evolutionary Computation, volume 23, number 2, pages 279-307, 2015.
> DOI: [10.1162/EVCO_a_00138](http://dx.doi.org/10.1162/EVCO_a_00138)

Before to use this package, code, and any resource attached to which,
and/or its derivatives, you must accept the **[license](LICENSE.md)**.
And, please, do not forget to refer the above paper.

[Contribution guidelines for this project](CONTRIBUTING.md)


Dependencies
-------------------------------

- [Boost libraries >= 1.58;](http://www.boost.org)

- Modern C++ compiler supporting C++11;
    - This code was tested using GCC 5.2 and GCC 7.2 successfully. Clang 3.9
      also compiles the code without problems, although changes in the
      parameters are necessary. A caveat is that to the date of this document
      was written, Clang does not use multiple threads of OpenMP, the
      technology behind the paralyzation of the decoding process.

- [IBM ILOG CPLEX >=
  12.5.0](https://www.ibm.com/products/ilog-cplex-optimization-studio) for
  algorithm versions depending on LP relaxations (CALP, GALP). CPLEX libs can
  be disabled on the Makefile;

This package also includes a modified version of the
[brkgaAPI](https://github.com/rfrancotoso/brkgaAPI). 
The core functionality was kept, but changes in the API make the distributed
version and the original version incompatible.

This package also includes (together brkgaAPI) a Mersenne Twister random number
generator licensed under BSD-3.


Compiling and testing
-------------------------------

First, you must set the proper paths to the "includes" and "libs" in your
system. In the `Makefile`, go to the Section "Lib and include definitions" 
(line 84) and

- In "Cplex and Concert settings" subsection (line 93), set the path to 
  your CPLEX installation. 
  
Note that depending of the configuration of your system, other adjustments
may be necessary.

To build the BRKGA, just use:
```bash
$ make all
```

