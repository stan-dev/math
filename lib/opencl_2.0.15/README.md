# OpenCL<sup>TM</sup> API C++ bindings

Doxgen documentation for the bindings is available here:

  http://khronosgroup.github.io/OpenCL-CLHPP/

Components:

  * `include/CL/opencl.hpp`:
    The latest, maintained, version of the C++ bindings. It should work with all
    versions of OpenCL (including 1.x). This is what most users will want.

  * `include/CL/cl2.hpp`:
    Includes `opencl.hpp` and emits a warning, for backwards compability.

  * `docs`:
    Doxygen file used to generate HTML documentation for `opencl.hpp`.

  * `examples`:
    A simple example application using the very basic features of the bindings.

  * `tests`:
    A (very small, incomplete) set of regression tests. Building the tests
    requires Python, Ruby, Unity and CMock. For the last two we use
    [Unity 2.1.0](https://github.com/ThrowTheSwitch/Unity/releases/tag/v2.1.0)
    and [CMock top-of-tree from Github](https://github.com/ThrowTheSwitch/CMock)
    (the version 2.0.204 on Sourceforge does not work).

  * `CMakeLists.txt`:
    Build system for the examples and tests and logic for the bindings
    installation.

To get external dependencies needed for testing, use `--recursive` when cloning
the repository, or run `git submodule update --init`.

You may need to tell CMake where to find the OpenCL headers and libraries,
using the variables `OPENCL_INCLUDE_DIR` and `OPENCL_LIB_DIR`.

These can be set either as environment variables, or on the cmake command line
using the syntax `-D<VAR>=<VALUE>`.

The following is an example set of commands to checkout and build the C++
bindings (adapt paths as required):

```
    git clone --recursive https://github.com/KhronosGroup/OpenCL-CLHPP
    cd OpenCL-CLHPP
    mkdir build
    cd build
    cmake .. -DOPENCL_INCLUDE_DIR=/path/to/OpenCL/headers -DOPENCL_LIB_DIR=/path/to/OpenCL/library
    make
    make test
```

After building, the headers appear in `build/include/CL/`.

If Doxygen is available, you can generate HTML documentation by typing `make docs`.
