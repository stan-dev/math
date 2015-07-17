The <b>Stan Math Library</b> is a C++, reverse-mode automatic differentiation library designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

Licensing
---------
The Stan Math Library is licensed under the [new BSD license](https://raw.githubusercontent.com/stan-dev/math/develop/licenses/math-license.txt).

Required Libraries
------------------
In order to run Stan Math, you *must* make available the following two libraries on which it depends:

- Boost (version 1.58): [Boost Home Page](http://www.boost.org)
- Eigen (version 3.24): [Eigen Home Page](http://eigen.tuxfamily.org/index.php?title=Main_Page)

These should be downloaded and unpacked into directories where they can be read.  

Only these two versions of the dependent libraries have been tested with Stan Math.

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
> clang++ -I /path/to/stan -I /path/to/Eigen -I /path/to/boost foo.cpp
> ./a.out
log normal(1 | 2, 3)=-2.07311
```

The `-I` includes point to the three necessary includes:

* the Stan Math library:  the path is to the source directory that contains `stan` as a subdirectory
* the Eigen C++ matrix library:  this path is to the source directory that contains `Eigen` as a subdirectory
* the Boost C++ library: this path is to the source directory that contains `boost` as a subdirectory

Note that the paths should *not* include the final directories `stan`, `Eigen`, or `boost` on the paths.  An example of a real instantiation:

```
clang++ -I ~/stan-dev/math -I ~/stan/lib/eigen_3.2.4/ -I ~/stan/lib/boost_1.58.0/ foo.cpp
```

The following directories all exist below the links given to `-I`: `~/stan-dev/math/stan` and `~/stan/lib/eigen_3.2.4/Eigen` and `~stan/lib/boost_1.58.0/boost`.

Other Compilers
---------------
There's nothing special about `clang++` --- the `g++` compiler behaves the same way.  You'll need to modify the commands for other compilers, which will need to be up-to-date enough to compile the Stan Math Library.
