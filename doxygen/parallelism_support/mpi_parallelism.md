# MPI {#mpi}

The [message passing interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) allows the exchange of messages between different processes. We can use MPI to parallelize the computation of a single log probability computation by using multiple processes. Since each process is *independent* of all other processes, we can run automatic differentiation (AD) in parallel by providing each process its own AD stack, then combining in a synchronized step which fits into the Math library's AD design.

The target audience for MPI are those with large computer clusters. For users looking for parallel computation on a single computer, please turn to a [[threading based approach|Threading Support]], which is easier to use and provides similar performance gains.

# Requirements {#mpi-requirements}

A base MPI installation must be installed on the system. See the [instructions from `boost.mpi`](http://www.boost.org/doc/libs/1_66_0/doc/html/mpi/getting_started.html#mpi.mpi_impl) to verify that there is a working MPI system.

Stan supports MPI on Mac OS X and Linux. **Windows is not supported at the moment.**

For Mac OS X and Linux, any MPI installation that works with `boost.mpi` is supported. The two major open-source base MPI implementations are `mpich` and `openMPI`. The Math library is tested with these two implementations while others supported by `boost.mpi` may work as well.

The base MPI installation provides the command line tools
+ `mpicxx`: The recommended  compiler command to use when building any MPI application.
+ `mpirun`: Wrapper binary used to start a MPI enabled binary on a given machine.

## Installation on a cluster system

Please ask your system administrator for details on how to compile, execute, and submit MPI applications.

## Installation on Mac OS X

Install `mpich` from MacPorts or homebrew.

## Installation on Linux

The package distribution system on your version of linux should have pre-built `mpich` (or `openmpi`) available.

In addition to that, you must have the following packages installed (Ubuntu package names listed): `python-dev, libxml2-dev, libxlst-dev`, and you may be required to add the following to your `make/local`: `LDLIBS+=-lpthread`.

## Installation on Windows

MPI is not supported on Windows at this time.

## Note on Boost

Stan builds it's own `boost.mpi` and `boost.serialization` libraries and installs these into its library subfolder. If the operating system provides these Boost libraries and it's required to use them, there is additional configuration that needs to be done (through `make/local`) to use that installation.

Moreover, the boost libraries are build using the boost build system. Boost build will attempt to auto-detect the MPI installation specifics on your system and the toolset to use. Should boost's auto-detect fail or a specific configuration be required, then users can configure the boost build system through the configuration file `stan-math/lib/boost_1.xx.x/user_config.jam` manually as needed.

## Note on compilers used

We strongly recommend to use the `mpicxx` command to build any program using MPI within Math. While it is possible to change the compiler used with these commands (openMPI has a `-cxx=` option, for example), this can only be done with great caution. The complication is that during compilation of the base MPI libraries the exact bit representation of each type is analyzed and strong deviations due to compiler changes may lead to unexpected behavior. In case of compiler mismatch between the base MPI libraries and `boost.mpi` (and Math) changes in the compiler ABI can lead to unexpected segfaults. Therefore, we recommend to use the `mpicxx` as compiler and do not recommend to deviate from the compiler used to build MPI. Often this means to use the system default compiler which may be rather old and not ideal for Stan. In such cases a more modern gcc (if gcc is the system compiler) can be considered as long as no ABI changes are known.

# Setting up the Math library with MPI

Stan uses the `boost.mpi` library to interface with the installed MPI implementation. `boost.mpi` is built automatically by the Math library when the Math library is configured for MPI. To configure MPI for the Math library, please proceed ass follows:

0. Ensure that a base MPI installation is available and accessible on the system. See [Requirements](#mpi-requirements).
1. Open a text file called `make/local`; if it does not exist, create one.
2. Add these lines to the `make/local` file:
```
STAN_MPI=true
CXX=mpicxx
```
3. Optional: instead of using `CXX=mpicxx`, the user can specify the compiler with the proper compiler and linker options needed to build an MPI enabled binary (the command `mpicxx -show` displays for `mpich` what is executed / `openmpi` uses `mpicxx -show-me`), but please read the note on compilers above.
4. Clean all binaries. After changing configuration through `make/local`, all the tests should be rebuilt. Please type:
```
make clean-all
```

Once the Math library is configured for MPI, the tests will be built with MPI. Note that the `boost.mpi` and `boost.serialization` library are build and linked against dynamically.

# Running tests with MPI

Once MPI is enabled, the `runTests.py` script in the `cmdstan/stan/lib/stan_math` directory will run all tests in an environment which resembles a MPI run. There are two types of tests:

- conventional tests: This includes all unit tests which do not use any MPI parallelism. In order to run these tests in a MPI like way we compile these with the `mpicxx` command and execute them with the `mpirun` run command. However, we explicitly disable for these serial tests the use of multiple processes. That is, the `runTests.py` script executes the tests with `mpirun -np 1 test/unit/.../.../some_test`. Starting up multiple threads for serial only tests would lead to race conditions since these codes are not prepared for parallelism.

- dedicated MPI tests: All tests matching the regular expression `*mpi_*test.cpp` will be executed by `runTests.py` with \`mpirun -np \#CPU test/unit/.../.../some_mpi_test.cpp\` and \`\#CPU\` will be set to the same argument as given to the `-j` option of `runTests.py`, but it will use at least 2 processes. This is to ensure that the MPI tests are actually run with multiple processes in parallel to emulate behavior under MPI execution. Note that `mpirun` is usually configured to disallow \#CPU to exceed the number of physical CPUs found on the machine.

To illustrate what is happening let's consider two examples (assuming MPI is enabled as described above):

- conventional test:
```
./runTests.py test/unit/math/prim/functor/map_rect_test.cpp
# => compilation with mpicxx
# => execution with mpirun using a single process
# mpirun -np 1 test/unit/math/prim/functor/map_rect_test
```

- dedicated MPI test:
```
./runTests.py test/unit/math/prim/functor/mpi_cluster_test.cpp
# => compilation with mpicxx
# => execution with mpirun using at least two processes
# mpirun -np 2 test/unit/math/prim/functor/mpi_cluster_test

./runTests.py -j8 test/unit/math/prim/functor/mpi_cluster_test.cpp
# => compilation with mpicxx
# => execution with mpirun using 8 processes
# mpirun -np 8 test/unit/math/prim/functor/mpi_cluster_test
```
