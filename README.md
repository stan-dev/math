The <b>Stan Math Library</b> is a C++, reverse-mode automatic differentiation library designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

Licensing
---------
The Stan Math Library is licensed under the [new BSD license](https://raw.githubusercontent.com/stan-dev/math/develop/licenses/stan-math-library-license.txt).

Required Libraries
------------------
Stan Math depends on two libraries:

- Boost (version 1.62.0): [Boost Home Page](http://www.boost.org)
- Eigen (version 3.2.9): [Eigen Home Page](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- Cvodes (version 2.9.0): [Sundials Home Page](http://computation.llnl.gov/projects/sundials/sundials-software)

These are distributed under the `lib/` subdirectory. Only these three versions of the dependent libraries have been tested with Stan Math.

Installation
------------
The Stan Math Library is a header-only C++ library.

A simple hello world program using Stan Math is as follows:

```
#include <stan/math.hpp>
#include <iostream>

int main() {
  std::cout << "log normal(1 | 2, 3)=" 
            << stan::math::normal_log(1, 2, 3) 
            << std::endl;
}
```

If this is in the file `/path/to/foo/foo.cpp`, then you can compile and run this with something like this, with the `path/to` business replaced with actual paths:

```
> cd /path/to/foo
> clang++ -I /path/to/stan-math -I /path/to/Eigen -I /path/to/boost foo.cpp  -I /path/to/cvodes
> ./a.out
log normal(1 | 2, 3)=-2.07311
```

The `-I` includes provide paths pointing to the four necessary includes:

* Stan Math Library:  path to source directory that contains `stan` as a subdirectory
* Eigen C++ Matrix Library:  path to source directory that contains `Eigen` as a subdirectory
* Boost C++ Library:  path to source directory that contains `boost` as a subdirectory
* Cvodes: path to source directory that contains `cvodes` as a subdirectory

Note that the paths should *not* include the final directories `stan`, `Eigen`, or `boost` on the paths.  An example of a real instantiation:

```
clang++ -I ~/stan-dev/math -I ~/stan-dev/math/lib/eigen_3.2.9/ -I ~/stan-dev/math/lib/boost_1.62.0/ -I ~/stan-dev/math/lib/cvodes_2.9.0/include foo.cpp
```

The following directories all exist below the links given to `-I`: `~/stan-dev/math/stan` and `~/stan-dev/math/lib/eigen_3.2.9/Eigen` and `~stan-dev/math/lib/boost_1.62.0/boost` and `~stan-dev/math/lib/cvodes_2.9.0/include`.

Other Compilers
---------------
There's nothing special about `clang++` --- the `g++` compiler behaves the same way.  You'll need to modify the commands for other compilers, which will need to be up-to-date enough to compile the Stan Math Library.
