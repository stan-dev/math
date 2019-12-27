The <b>Stan Math Library</b> is a C++, reverse-mode automatic differentiation library designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

[![DOI](https://zenodo.org/badge/38388440.svg)](https://zenodo.org/badge/latestdoi/38388440)

Licensing
---------
The Stan Math Library is licensed under the [new BSD license](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd).

The Stan Math Library depends on the Intel TBB library which is licensed under the Apache 2.0 license. This dependency implies an additional restriction as compared to the new BSD lincense alone. The Apache 2.0 license is incompatible with GPL-2 licensed code if distributed as a unitary binary. You may refer to the Apache 2.0 evaluation page on the [Stan Math wiki](https://github.com/stan-dev/math/wiki/Apache-2.0-License-Evaluation).

Required Libraries
------------------
Stan Math depends on four libraries:

- Boost (version 1.69.0): [Boost Home Page](https://www.boost.org)
- Eigen (version 3.3.3): [Eigen Home Page](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- SUNDIALS (version 4.1.0): [Sundials Home Page](https://computation.llnl.gov/projects/sundials/sundials-software)
- Intel TBB (version 2019_U8): [Intel TBB Home Page](https://www.threadingbuildingblocks.org)

These are distributed under the `lib/` subdirectory. Only these versions of the dependent libraries have been tested with Stan Math.

Documentation
------------

Documentation for Stan math is available at [mc-stan.org/math](https://mc-stan.org/math/)

Installation
------------
The Stan Math Library is a C++ library which depends on the Intel TBB library and requires for some functionality (ordinary differential equations and root solving) on the Sundials library. As build system the make facility is used to manage all dependencies.

A simple hello world program using Stan Math is as follows:

```cpp
#include <stan/math.hpp>
#include <iostream>

int main() {
  std::cout << "log normal(1 | 2, 3)="
            << stan::math::normal_log(1, 2, 3)
            << std::endl;
}
```

If this is in the file `/path/to/foo/foo.cpp`, then you can compile and run this with something like this, with the `path/to` business replaced with actual paths:

```bash
> cd /path/to/foo
> make -f path/to/stan-math/make/standalone math-libs
> make -f path/to/stan-math/make/standalone foo
> ./foo
log normal(1 | 2, 3)=-2.07311
```

The first make command with the `math-libs` target ensures that all binary dependencies are build and ready to use. The second make command ensures that all dependencies of Stan Math are available to the compiler which are:

* Stan Math Library:  path to source directory that contains `stan` as a subdirectory
* Eigen C++ Matrix Library:  path to source directory that contains `Eigen` as a subdirectory
* Boost C++ Library:  path to source directory that contains `boost` as a subdirectory
* SUNDIALS: path to source directory that contains `cvodes` and `idas` as a subdirectory
* Intel TBB: path to source directory that contains `tbb` as a subdirectory

Note that the paths should *not* include the final directories `stan`, `Eigen`, or `boost` on the paths.  An example of a real instantiation:

```bash
> make -f ~/stan-dev/math/make/standalone math-libs
> make -f ~/stan-dev/math/make/standalone foo
```

The `math-libs` target must only be called a single time and can be omitted for subsequent compilations.

The standalone makefile ensure that all the required `-I` include statements are given to the compiler and the necessary libraries are linked: `~/stan-dev/math/stan` and `~/stan-dev/math/lib/eigen_3.3.3/Eigen` and `~/stan-dev/math/lib/boost_1.69.0/boost` and `~/stan-dev/math/lib/sundials_4.1.0/include` and `~/stan-dev/math/lib/tbb_2019_U8/include`. The `~/stan-dev/math/lib/tbb` directory is created by the math-libs makefile target automatically. The `-Wl,-rpath,...` instruct the linker to hard-code the path to the Intel TBB library inside the stan-math directory into the final binary. This way the Intel TBB is found when executing the program.

Note for Windows users: On Windows the rpath feature as used by Stan Math to hardcode an absolute path to a dynamically loaded library does not work. On Windows the TBB dynamic library `tbb.dll` is located in the `math/lib/tbb` directory. The user can choose to either copy this file to the same directory of the executable or to add the directory `path/to/math/lib/tbb` as absolute path to the system wide `PATH` variable.

Other Compilers
---------------
There's nothing special about `clang++` --- the `g++` compiler behaves the same way.  You'll need to modify the commands for other compilers, which will need to be up-to-date enough to compile the Stan Math Library.
