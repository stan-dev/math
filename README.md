<div><a href="https://zenodo.org/badge/latestdoi/38388440"><img src="https://zenodo.org/badge/38388440.svg"/></a></div>

The <b>Stan Math Library</b> is a C++, reverse-mode automatic
differentiation library designed to be usable, extensive and
extensible, efficient, scalable, stable, portable, and redistributable
in order to facilitate the construction and utilization of algorithms
that utilize derivatives.


Licensing
---------
The Stan Math Library is licensed under the [new BSD
license](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd).

The Stan Math Library depends on the Intel TBB library which is
licensed under the Apache 2.0 license. This dependency implies an
additional restriction as compared to the new BSD license alone. The
Apache 2.0 license is incompatible with GPL-2 licensed code if
distributed as a unitary binary. You may refer to the Licensing page on the [Stan wiki](https://github.com/stan-dev/stan/wiki/Stan-Licensing).

Required Libraries
------------------
Stan Math depends on four libraries:

- Boost (version 1.78.0): [Boost Home Page](https://www.boost.org)
- Eigen (version 3.3.9: [Eigen Home Page](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- SUNDIALS (version 6.1.1): [Sundials Home Page](https://computing.llnl.gov/projects/sundials)
- Intel TBB (version 2020.3): [Intel TBB Home Page](https://www.threadingbuildingblocks.org)

These are distributed under the `lib/` subdirectory. Only these
versions of the dependent libraries have been tested with Stan Math.

Documentation
------------

Documentation for Stan math is available at [mc-stan.org/math](https://mc-stan.org/math/)

Contributing
------------

We love contributions from everyone in the form of good discussion, issues, and pull requests.
If you are interested in contributing to Stan math please check the Contributor Guide at [mc-stan.org/math](https://mc-stan.org/math/).

Installation
------------
The Stan Math Library is a C++ library which depends on the Intel TBB
library and requires for some functionality (ordinary differential
equations and root solving) the Sundials library. The build system is
the make facility, which is used to manage all dependencies.

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

If this is in the file `/path/to/foo/foo.cpp`, then you can compile
and run this with something like this, with the `/path/to` business
replaced with actual paths:

```bash
> cd /path/to/foo
> make -j4 -f /path/to/stan-math/make/standalone math-libs
> make -f /path/to/stan-math/make/standalone foo
> ./foo
log normal(1 | 2, 3)=-2.07311
```

The first make command with the `math-libs` target ensures that all
binary dependencies of Stan Math are built and ready to use. The `-j4`
instructs `make` to use 4 cores concurrently which should be adapted
to your needs. The second make command ensures that the Stan Math
sources and all of the dependencies are available to the compiler when
building `foo`.

An example of a real instantiation whenever the path to Stan Math is
`~/stan-dev/math/`:

```bash
> make -j4 -f ~/stan-dev/math/make/standalone math-libs
> make -f ~/stan-dev/math/make/standalone foo
```
The `math-libs` target has to be called only once, and can be omitted
for subsequent compilations.

The standalone makefile ensures that all the required `-I` include
statements are given to the compiler and the necessary libraries are
linked: `~/stan-dev/math` and `~/stan-dev/math/lib/eigen_3.3.9` and
`~/stan-dev/math/lib/boost_1.78.0` and
`~/stan-dev/math/lib/sundials_6.1.1/include` and
`~/stan-dev/math/lib/tbb_2020.3/include`. The
`~/stan-dev/math/lib/tbb` directory is created by the `math-libs`
makefile target automatically and contains the dynamically loaded
Intel TBB library. The flags `-Wl,-rpath,...` instruct the linker to
hard-code the path to the dynamically loaded Intel TBB library inside
the stan-math directory into the final binary. This way the Intel TBB
is found when executing the program.

Note for Windows users: On Windows the `-rpath` feature as used by
Stan Math to hardcode an absolute path to a dynamically loaded library
does not work. On Windows the Intel TBB dynamic library `tbb.dll` is
located in the `math/lib/tbb` directory. The user can choose to copy
this file to the same directory of the executable or to add the
directory `/path/to/math/lib/tbb` as absolute path to the system-wide
`PATH` variable.

Intel TBB
---------

`math` supports the new interface of Intel TBB, can be configured to use an external copy of TBB (e.g., with [`oneTBB`](https://github.com/oneapi-src/oneTBB) or the system TBB library), using the `TBB_LIB` and `TBB_INC` environment variables.

To build the development version of `math` with [`oneTBB`](https://github.com/oneapi-src/oneTBB):

- Install [`oneTBB`](https://github.com/oneapi-src/oneTBB).

For example, installing [`oneTBB`](https://github.com/oneapi-src/oneTBB) on Linux 64-bit (`x86_64`) to `$HOME` directory (change if needed!):
```bash
TBB_RELEASE="https://api.github.com/repos/oneapi-src/oneTBB/releases/latest"
TBB_TAG=$(curl --silent $TBB_RELEASE | grep -Po '"tag_name": "\K.*?(?=")')
TBB_VERSION=${TBB_TAG#?}

wget https://github.com/oneapi-src/oneTBB/releases/download/v${TBB_VERSION}/oneapi-tbb-${TBB_VERSION}-lin.tgz
tar zxvf oneapi-tbb-$TBB_VERSION-lin.tgz -C $HOME

export TBB="$HOME/oneapi-tbb-$TBB_VERSION"
```
Note that you may replace `TBB_VERSION=${TBB_TAG#?}` with a custom version number if needed ( check available releases [here](https://github.com/oneapi-src/oneTBB/releases) ).

- Set the TBB environment variables (specifically: `TBB` for the installation prefix, `TBB_INC` for the directory that includes the header files, and `TBB_LIB` for the libraries directory).

For example, installing [`oneTBB`](https://github.com/oneapi-src/oneTBB) on Linux 64-bit (`x86_64`) to `$HOME` directory (change if needed!):
```bash
source $TBB/env/vars.sh intel64

export TBB_INC="$TBB/include"
export TBB_LIB="$TBB/lib/intel64/gcc4.8"
```

- Set `Stan` local compiler flags to use the new TBB interface:
```bash
mkdir -p ~/.config/stan
echo TBB_INTERFACE_NEW=true>> ~/.config/stan/make.local
```

Compilers
---------

The above example will use the default compiler of the system as
determined by `make`. On Linux this is usually `g++`, on MacOS
`clang++`, and for Windows this is `g++` if the RTools for Windows are
used. There's nothing special about any of these and they can be
changed through the `CXX` variable of `make`. The recommended way to
set this variable for the Stan Math library is by creating a
`make/local` file within the Stan Math library directory. Defining
`CXX=g++` in this file will ensure that the GNU C++ compiler is always
used, for example. The compiler must be able to fully support C++11
and partially the C++14 standard. The `g++` 4.9.3 version part of
RTools for Windows currently defines the minimal C++ feature set
required by the Stan Math library.

Note that whenever the compiler is changed, the user usually must
clean and rebuild all binary dependencies with the commands:
```bash
> make -f path/to/stan-math/make/standalone math-clean
> make -j4 -f path/to/stan-math/make/standalone math-libs
```
This ensures that the binary dependencies are created with the new
compiler.
