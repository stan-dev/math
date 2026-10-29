<div><a href="https://zenodo.org/badge/latestdoi/38388440"><img src="https://zenodo.org/badge/38388440.svg"/></a></div>

The **Stan Math Library** is a C++ library for reverse-mode automatic
differentiation. Use Stan Math to build algorithms requiring derivatives.


Documentation, Installation, and Examples
--------------

Read the Stan Math documentation at mc-stan.org.
The wiki articles are outdated. Stan Math developers do not maintain wiki articles.

[mc-stan.org/math](https://mc-stan.org/math/)


Licensing
---------
The Stan Math Library is licensed under the [new BSD
license](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd).

Stan Math depends on Intel TBB. TBB uses the Apache 2.0 license.
The Apache 2.0 license adds restrictions to the BSD license terms.
The Apache 2.0 license is incompatible with GPL-2 licensed code in a combined binary. See the [Licensing page](https://github.com/stan-dev/stan/wiki/Stan-Licensing) for details.

Building tests with CMake
-------------------------

Stan Math has recently used to using cmake as its primary build system.
CMake creates build rules for a selected build tool while CTest is CMake's test runner.

`runTests.py` is the main way to compile the tests for Stan Math and it now uses cmake on the backend.
`runTests.py` configures CMake, builds the selected executables, and runs the tests.

```sh
# Builds a single executable for the whole `fun` directory
./runTests.py test/unit/math/prim/fun
# Builds a single executable for the `sum_test.cpp` file
./runTests.py test/unit/math/prim/fun/sum_test.cpp
```

Some useful options for `runTests.py` include:

- `--fresh` resets cached CMake settings before configuration.

CMake saves persistent build settings in a cache.
The default cache file is `build/CMakeCache.txt`.
Cached settings can override later changes to `make/local`.
`--fresh` resets the CMake configuration before building.
`--fresh` removes `CMakeCache.txt` and the top-level `CMakeFiles` directory, which is what CMake's own `--fresh` option does.
It works with any supported CMake version, while the native option needs CMake 3.24 or newer.
CMake then configures the build with the current arguments and `make/local` settings.
Objects whose compile flags did not change are not rebuilt.
Use `--fresh` when cached settings no longer match the intended build.
Pass required custom settings again because `--fresh` removes cached overrides.
`--test-only` runs an existing executable without configuration or compilation.
`--fresh` cannot be combined with `--test-only`.

```sh
# Reset cached settings, then build and run the sum tests.
./runTests.py --fresh test/unit/math/prim/fun/sum_test.cpp
```

- `--build-dir PATH` selects the build directory. The default is `build`.
The build directory contains generated build rules, compiled files, and test results.

```sh
# Store the build in `.build`.
./runTests.py --build-dir .build test/unit/math/prim/fun/sum_test.cpp
```

- `--config NAME` selects the build configuration such as `Debug`, `Release`, or `ReleaseTesting`. Setting a configuration changes the default behavior of the build system such as compiler flags and debug information.
The default is `ReleaseTesting`.
`ReleaseTesting` compiles with `-O1` and does not define `NDEBUG`, so assertions stay on as in the make build.
`Release` defines `NDEBUG`, which disables assertions.
`STAN_TEST_OPTIMIZATION` still sets the optimization level of the tests in every configuration.

```sh
# Use the `Release` configuration.
./runTests.py --config Release test/unit/math/prim/fun/sum_test.cpp
```

- `-j N` sets the number of parallel build jobs. The default is 1.

```sh
# Build with four jobs.
./runTests.py -j4 test/unit/math/prim/fun
```

- `-f TEXT` selects source files whose paths contain `TEXT`.
A directory executable includes only source files matching the filter.
Repeat `-f` to require all filters to match.

```sh
# Build and run primitive function tests with `sum` in the source paths.
./runTests.py -f sum test/unit/math/prim/fun
```

- `--make-only` builds the tests without running the test executables.

```sh
# Build the executable for `sum_test.cpp`.
./runTests.py --make-only test/unit/math/prim/fun/sum_test.cpp
```

Selected files and distribution fixtures build through one CMake target, `stan_selected_tests`, so all of them compile in parallel.
Each selected file is still its own executable.

- A distribution fixture such as `test/prob/normal/normal_test.hpp` selects the tests generated from that fixture.

```sh
# Build and run the tests generated from the normal distribution fixture.
./runTests.py -j8 test/prob/normal/normal_test.hpp
```

- `--configure-only` configures `<build-dir>` without building or running tests.
`test/prob/getDependencies.py` reads that configuration's compile commands to select the distribution fixtures whose included Stan headers changed:

```sh
./runTests.py --configure-only
./runTests.py -j8 $(python3 test/prob/getDependencies.py --build-dir build)
```

- `--test-only` runs an existing test executable.
`--test-only` skips configuration and compilation.
`--test-only` does not rebuild tests after source changes.

```sh
# Build the primitive function test executable in `.build`
./runTests.py --make-only --build-dir .build test/unit/math/prim/fun
# run the executable.
./runTests.py --test-only --build-dir .build test/unit/math/prim/fun
```

- `--changed` selects tests for committed changes on your branch.
`--changed` compares your branch with `origin/develop`.
`--changed` selects changed test files and tests for changed Stan headers.
Each selected file builds and runs as a separate executable:

```sh
./runTests.py --changed
```

- `--run-all` continues to the next executable if a test executable fails.
`runTests.py --run-all` still returns a failure status if any tests fail.
Build failures stop `runTests.py`.

```sh
# Run the `abs` tests even if the `sum` tests fail.
./runTests.py --run-all test/unit/math/prim/fun/sum_test.cpp test/unit/math/prim/fun/abs_test.cpp
```

- `--test-jobs N` runs up to `N` test executables at once. The default is 1.
Each executable prints its output when it finishes, in the order of the selected paths.
With more than one job, every executable runs even if one fails.
A trailing `/` in a pattern selects only directories.

```sh
# Build one executable per distribution directory and run four at a time.
./runTests.py -j8 --test-jobs 4 'test/prob/*/'
```

`<build-dir>` means the directory selected by `--build-dir`.
XML results are written under `<build-dir>/test-results` by default.
`make/local` is an optional file for local compiler and build settings.
The CMake test setup reads `make/local`.
All CMake generators use the imported settings, including Ninja and Unix Makefiles.
The Windows and ARM GitHub jobs use CMake.

`--cmake-arg=ARGUMENT` passes one argument to CMake configuration.
Repeat `--cmake-arg` to pass multiple arguments.
The CMake argument `-DNAME=VALUE` sets a named cache setting.

```sh
# Build reverse-mode core tests with threading enabled.
./runTests.py --cmake-arg=-DSTAN_THREADS=ON test/unit/math/rev/core
```

`-DSTAN_SANITIZE_ADDRESS=ON` compiles C++ sources and links with AddressSanitizer.
The vendored TBB build does not use it.

```sh
# Build mixed-mode tests at -O1 with AddressSanitizer, without make/local.
./runTests.py --cmake-arg=-DSTAN_USE_MAKE_LOCAL=OFF --cmake-arg=-DSTAN_SANITIZE_ADDRESS=ON \
  --cmake-arg=-DSTAN_TEST_OPTIMIZATION=1 test/unit/math/mix
```

The `cpplint` and `test-math-dependencies` CMake targets run the same checks as `make cpplint` and `make test-math-dependencies`:

```sh
cmake -S . -B build -DBUILD_TESTING=OFF
cmake --build build --target cpplint
```

The CMake argument `-G NAME` selects a generator.
The combined form `-GNAME` selects the generator with one argument.

```sh
# Build with Ninja and the settings from `make/local`.
./runTests.py --build-dir build-ninja --cmake-arg=-GNinja test/unit/math/prim/fun/sum_test.cpp
```

CMake imports the following Make settings:

| Make setting | CMake setting |
| --- | --- |
| `CC`, `CXX` | C and C++ compilers |
| `O` | Optimization level from 0 to 3, mapped to `STAN_TEST_OPTIMIZATION` |
| `STAN_THREADS` | Enable threading support |
| `STAN_OPENCL` | Enable OpenCL GPU support |
| `STAN_MPI` | Enable MPI support for communication between processes |
| `STAN_NO_RANGE_CHECKS` | Disable range checks |
| `STAN_TEST_ROW_VECTORS` | Include row vectors in distribution tests |
| `STAN_PROB_TEST_ALL` | Enable all distribution test modes |
| `OPENCL_DEVICE_ID` | Select an OpenCL device with `STAN_OPENCL_DEVICE_ID` |
| `OPENCL_PLATFORM_ID` | Select an OpenCL platform with `STAN_OPENCL_PLATFORM_ID` |
| `INTEGRATED_OPENCL` | Enable integrated GPU settings with `STAN_INTEGRATED_OPENCL` |
| `BOOST`, `EIGEN`, `TBB`, `SUNDIALS`, `OPENCL`, `BENCHMARK` | Dependency source directories |
| `CXXFLAGS`, `CXXFLAGS_OPTIM`, `CXXFLAGS_WARNINGS` | Additional C++ compiler flags |
| `CFLAGS` | Additional C compiler flags |
| `CPPFLAGS`, `CPPFLAGS_OPTIM` | Additional C and C++ preprocessor flags |
| `LDFLAGS`, `LDFLAGS_OPTIM` | Additional linker flags |
| `LDLIBS` | Additional link libraries |

CMake options passed with `-D` take priority over imported defaults.
Existing CMake cache values also take priority over imported defaults.
Use a new build directory after changing cached settings or compilers.
Use `--fresh` to discard cached settings in the selected build directory.
Flags from `make/local` are added to the CMake compiler and linker flags.
Changes to `make/local` and included Make files trigger CMake configuration.
Make settings outside the table are not imported.

`CXX = clang++` selects the Clang C++ compiler.
`O = 0` disables compiler optimization.
`STAN_THREADS = true` enables threading support.
`CXXFLAGS += -g` adds compiler debug information.
Example: Set compiler settings in `make/local`.

```make
CXX = clang++
O = 0
STAN_THREADS = true
CXXFLAGS += -g
```

Use `STAN_USE_MAKE_LOCAL=OFF` to ignore `make/local`:

```sh
./runTests.py --cmake-arg=-DSTAN_USE_MAKE_LOCAL=OFF test/unit/math/prim/fun/sum_test.cpp
```

### Common CMake commands

Run the commands from the repository root.
CMake configuration creates build rules and saves settings in `{BUILD_DIRECTORY}/CMakeCache.txt`.
Configuration does not build the test executables.
`cmake --build` compiles and links the test executables.
CTest runs test executables already built by CMake.
`runTests.py` configures, builds, and runs the selected tests automatically.

#### Configure a build

`-S PATH` selects the source directory. Use `.` for the repository root.
`-B PATH` selects the build directory. CMake creates the directory if needed.
`-DNAME=VALUE` sets a named cache setting.
`CMAKE_BUILD_TYPE` selects the configuration for Ninja and Unix Makefiles.
`Release` enables release compiler settings.
`ReleaseTesting` is the configuration `runTests.py` uses: `-O1` without `NDEBUG`.
Direct CMake commands also read `make/local`.
Example: Configure a release build in `build`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

The native CMake option `--fresh` resets the configuration cache.
`--fresh` requires CMake 3.24 or newer.
With older CMake, delete the cache file and the top-level `CMakeFiles` directory instead.
Example: Reset cached settings and configure a release build.

```sh
cmake --fresh -S . -B build -DCMAKE_BUILD_TYPE=Release
# CMake older than 3.24
cmake -E remove -f build/CMakeCache.txt
cmake -E remove_directory build/CMakeFiles
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

#### List settings and descriptions

`-N` reads the cache without configuring the project.
`-L` lists cache variables.
`H` in `-LH` or `-LAH` prints the description above each value.
`A` in `-LAH` includes advanced variables.
Advanced variables contain settings changed less often during normal builds.
Cache listings omit temporary variables and uncached `make/local` defaults.
Print the current cache values and option descriptions:

```sh
cmake -N -LH build
```

Add `A` to include advanced cache settings:

```sh
cmake -N -LAH build
```

To inspect project defaults, configure an empty build directory.
Disable `make/local` imports to exclude local overrides.
The `build-defaults` directory must not contain an existing CMake cache:

```sh
cmake -S . -B build-defaults -DSTAN_USE_MAKE_LOCAL=OFF -LAH
```

Project option definitions contain the default values and descriptions.
Read [CMakeLists.txt](CMakeLists.txt) for the main options.
Read [cmake/CMakeLists.txt](cmake/CMakeLists.txt) for dependency options.
Read [cmake/MakeLocal.cmake](cmake/MakeLocal.cmake) for local settings options.
Read [test/prob/CMakeLists.txt](test/prob/CMakeLists.txt) for distribution options.

`--help-variable-list` prints CMake's built-in variable names.
`--help-variable NAME` prints the definition of a built-in variable.
CMake's built-in help does not define Stan-specific options.
Use the cache descriptions or the project CMake files for Stan-specific options.
Example: List built-in variables and read the definition of `CMAKE_BUILD_TYPE`.

```sh
cmake --help-variable-list
cmake --help-variable CMAKE_BUILD_TYPE
```

#### Change settings or select a generator

Use `-DNAME=VALUE` to change a cached setting.
Example: Enable threading in `build`.

```sh
cmake -S . -B build -DSTAN_THREADS=ON
```

Use `-U` to remove a named cache override.
Configuration then uses the `make/local` value or the project default.
Example: Remove the cached threading setting.

```sh
cmake -S . -B build -U STAN_THREADS
```

Run `cmake --help` to list the available generators.
Use `-G` to select Ninja or Unix Makefiles.
Use a separate build directory when changing generators or compilers:

```sh
cmake --help
cmake -S . -B build-ninja -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake -S . -B build-make -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release
```

#### List targets and build tests

A build target names an executable, library, or build task.
`--build PATH` builds targets from the selected build directory.
`--target NAME` selects a build target.
The `help` target lists available targets for Ninja and Unix Makefiles.
`test_unit_math_prim_fun_sum_test` builds the executable for `sum_test.cpp`.
`--parallel N` sets the number of parallel build jobs.
`--verbose` prints compiler and linker commands.
List build targets for Ninja or Unix Makefiles:

```sh
cmake --build build --target help
```

Build one test executable with four parallel jobs:

```sh
cmake --build build --parallel 4 --target test_unit_math_prim_fun_sum_test
```

Add `--verbose` to print compiler and linker commands:

```sh
cmake --build build --parallel 4 --verbose --target test_unit_math_prim_fun_sum_test
```

#### List and run tests

CTest does not build missing executables.
Build the selected test targets before running CTest.
The CTest option `--test-dir PATH` selects the build directory containing test registrations.
A regular expression is a pattern used to match text.
The CTest option `-N` lists registered tests without running the tests.
`-R PATTERN` selects test names matching a regular expression.
`^` matches the start of a name. `$` matches the end of a name.
`--output-on-failure` prints output from failing tests.
The CTest option `--parallel N` sets the number of parallel test jobs.
List registered tests without running the test executables:

```sh
ctest --test-dir build -N
```

Run one registered test and print output on failure:

```sh
ctest --test-dir build -R '^test_unit_math_prim_fun_sum_test$' --output-on-failure
```

Run all registered tests with four parallel jobs:

```sh
ctest --test-dir build --parallel 4 --output-on-failure
```

See the [CMake command reference](https://cmake.org/cmake/help/latest/manual/cmake.1.html)
and [CTest command reference](https://cmake.org/cmake/help/latest/manual/ctest.1.html)
for more commands.

Each test source has a separate executable and CTest entry.
Targets including `unit_math_rev_core_subtests` and `prob_math_bernoulli_subtests`
build tests in the selected directory and the selected directory's subdirectories.
Use `-R` to select matching CTest entries:

```sh
ctest --test-dir build -R '^test_unit_math_rev_core_' --output-on-failure
```

The optional `*_all_tests` targets combine a directory's tests into one executable.
The default build and CTest exclude `*_all_tests` targets to prevent duplicate test runs.
Each test source compiles into one object file.
Individual test executables and directory executables reuse the compiled object files.
Switching between individual and directory executables does not recompile the test source.
The working directory for tests is the repository root.

CMake creates targets for all enabled tests.
Select a directory target to build only the required test executables.
The `unit_math_prim_core_subtests` target builds primitive core test executables.
The CTest label `test/unit/math/prim/core/` selects primitive core tests.
The CTest option `-L PATTERN` selects labels matching a regular expression.

```sh
cmake --build build --parallel 4 --target unit_math_prim_core_subtests
ctest --test-dir build -L '^test/unit/math/prim/core/' --output-on-failure
```

Precompiled headers (PCH) are enabled by default with `STAN_TEST_PCH=ON`.
Upstream Clang 9 and 10 automatically compile without PCH because of a compiler
defect fixed in Clang 11. Configuration reports the requested and effective
settings; the same setting applies to var-matrix compatibility checks.
Automatic differentiation (AD) computes derivatives through the Stan Math type system.
Unit tests reuse a PCH for the automatic differentiation layer under test.
Sources including `test_ad.hpp` use the AD helper PCH instead.
The PCH choice for each source also applies to directory executables.
Sources changing Stan or Eigen configuration macros skip PCH.
Sources without configuration macro changes continue to use PCH.
Umbrella header checks and tests with multiple translation units compile without PCH.

Distribution tests precompile distribution fixture headers.
The distribution test PCH excludes generated typedefs.
PCH builds and test builds use identical compiler flags and definitions.
GCC reports an invalid PCH as an error.
Use `-DSTAN_TEST_PCH=OFF` to compile tests without PCH.
Standalone header targets never use PCH.
Examples include `test-headers` and `stan_math_prim_fun_sum-test`.

Tests and the included SUNDIALS library use `STAN_TEST_OPTIMIZATION`.
`STAN_TEST_OPTIMIZATION` defaults to `1`. Optimization level `1` matches Make builds with `O=1`.

The build generates distribution test sources under `<build-dir>/test/prob`.
The distribution test generator is `test/prob/generate_tests.cpp`.
`STAN_PROB_TEST_BATCH_SIZE` defaults to 100 cases per source.
The build generates distribution test sources once per fixture before compilation.
The build regenerates distribution test sources when a fixture or the generator changes.

`STAN_THREADS`, `STAN_OPENCL`, `STAN_MPI`, and `STAN_NO_RANGE_CHECKS` work as before.
The build excludes OpenCL tests when `STAN_OPENCL=OFF`.
The build excludes MPI tests when `STAN_MPI=OFF`.
Use `STAN_USE_SYSTEM_{BOOST,EIGEN,TBB,SUNDIALS}` to select installed dependencies.
Use `STAN_{BOOST,EIGEN,TBB,SUNDIALS,OPENCL,BENCHMARK}_ROOT` to set dependency paths.
The included TBB source supports GCC and Clang builds on Linux and macOS.
MinGW and Rtools builds on Windows use the included TBB source.
MSVC builds require an installed TBB.
MPI builds use the included Boost sources.
Set `STAN_USE_SYSTEM_BOOST=ON` to use installed Boost libraries.
Set `BUILD_TESTING=OFF` to use the `stanmath` interface target without configuring tests.

Jenkins uses CMake and CTest. Each worker reuses the `build` directory.
Jenkins passes `--fresh` when a configuration cache already exists.
Test reports are saved under `test-results` with a directory for each test stage.
`PARALLEL` sets the number of build jobs.
The pipeline includes ASan, GPU, MPI, threading, and expression tests.
`STAN_TEST_ROW_VECTORS` and `STAN_PROB_TEST_ALL` enable the additional
distribution cases requested by Jenkins.
`runTests.py generate-expressions PATH` creates expression test sources in
`PATH/expression-sources`.
`STAN_BUILD_EXPRESSION_TESTS=ON` enables generated expression tests.
`runTests.py` also generates expression sources when you select `test/expressions`.
The runner option `-e N` sets the number of expression source files.
The runner option `--only-functions` selects functions for expression tests.
Example: Generate expression test sources in `build-expressions/expression-sources`.

```sh
python3 runTests.py generate-expressions build-expressions
```

Expression framework self-tests use the CMake runner.
Var-matrix compatibility checks compile through CMake.
The Makefiles remain available for upstream packages.
