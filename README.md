The <b>Stan Math Library</b> is a C++, reverse-mode automatic differentiation library designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

[![DOI](https://zenodo.org/badge/38388440.svg)](https://zenodo.org/badge/latestdoi/38388440)

Licensing
---------
The Stan Math Library is licensed under the [new BSD license](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd).

The Stan Math Library depends on the Intel TBB library which is licensed under the Apache 2.0 license. This dependency implies an additional restriction as compared to the new BSD lincense alone. The Apache 2.0 license is incompatible with GPL-2 licensed code if distributed as a unitary binary. You may refer to the Apache 2.0 evaluation page on the [Stan Math wiki](https://github.com/stan-dev/math/wiki/Apache-2.0-License-Evaluation).

Required Libraries
------------------
Stan Math depends on four libraries:

- Boost (version 1.72.0): [Boost Home Page](https://www.boost.org)
- Eigen (version 3.3.3): [Eigen Home Page](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- SUNDIALS (version 4.1.0): [Sundials Home Page](https://computation.llnl.gov/projects/sundials/sundials-software)
- Intel TBB (version 2019_U8): [Intel TBB Home Page](https://www.threadingbuildingblocks.org)

These are distributed under the `lib/` subdirectory. Only these versions of the dependent libraries have been tested with Stan Math.

Documentation
------------

Documentation for Stan math is available at [mc-stan.org/math](https://mc-stan.org/math/)

Installation
------------
The Stan Math Library is a C++ library which depends on the Intel TBB library and requires for some functionality (ordinary differential equations and root solving) on the Sundials library. The build system used is the make facility, which is used to manage all dependencies.

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

The first make command with the `math-libs` target ensures that all binary dependencies are built and ready to use. The second make command ensures that all dependencies of Stan Math are available to the compiler. These are:

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
The `math-libs` target has to be called only once, and can be omitted for subsequent compilations.

The standalone makefile ensures that all the required `-I` include statements are given to the compiler and the necessary libraries are linked: `~/stan-dev/math/stan` and `~/stan-dev/math/lib/eigen_3.3.3/Eigen` and `~/stan-dev/math/lib/boost_1.69.0/boost` and `~/stan-dev/math/lib/sundials_4.1.0/include` and `~/stan-dev/math/lib/tbb_2019_U8/include`. The `~/stan-dev/math/lib/tbb` directory is created by the `math-libs` makefile target automatically. The flags `-Wl,-rpath,...` instruct the linker to hard-code the path to the Intel TBB library inside the stan-math directory into the final binary. This way the Intel TBB is found when executing the program.

Note for Windows users: On Windows the `-rpath` feature as used by Stan Math to hardcode an absolute path to a dynamically loaded library does not work. On Windows the Intel TBB dynamic library `tbb.dll` is located in the `math/lib/tbb` directory. The user can choose to copy this file to the same directory of the executable or to add the directory `path/to/math/lib/tbb` as absolute path to the system-wide `PATH` variable.

Compilers
---------
The above example will use the default compiler of the system as determined by `make`. On Linux this is usually `g++`, on MacOS `clang++`, and for Windows this is `g++` if the RTools for Windows are used. There's nothing special about any of these and they can be changed through the `CXX` variable of `make`. The recommended way to set this variable for the Stan Math library is by creating a `make/local` file within the Stan Math library directory. Defining `CXX=g++` in this file will ensure that the GNU C++ compiler is always used, for example. The compiler must be able to fully support C++11 and partially the C++14 standard. The `g++` 4.9.3 version part of RTools for Windows currently defines the minimal C++ feature set required by the Stan Math library.

Note that whenever the compiler is changed, the user usually must clean and rebuild all binary dependencies with the commands:
```bash
> make -f path/to/stan-math/make/standalone math-clean
> make -f path/to/stan-math/make/standalone math-libs
```
This ensures that the binary dependencies are created with the new compiler.
