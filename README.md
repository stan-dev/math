<div><a href="https://zenodo.org/badge/latestdoi/38388440"><img src="https://zenodo.org/badge/38388440.svg"/></a></div>

The <b>Stan Math Library</b> is a C++, reverse-mode automatic
differentiation library designed to be usable, extensive and
extensible, efficient, scalable, stable, portable, and redistributable
in order to facilitate the construction and utilization of algorithms
that utilize derivatives.


Documentation, Installation, and Examples
--------------

All of Stan math's documentation is hosted on our website below. Please do not
reference articles in the wiki as they are outdated and not maintained.

[mc-stan.org/math](https://mc-stan.org/math/)


Licensing
---------
The Stan Math Library is licensed under the [new BSD
license](https://github.com/stan-dev/math/blob/develop/LICENSE%2Emd).

The Stan Math Library depends on the Intel TBB library which is
licensed under the Apache 2.0 license. This dependency implies an
additional restriction as compared to the new BSD license alone. The
Apache 2.0 license is incompatible with GPL-2 licensed code if
distributed as a unitary binary. You may refer to the Licensing page on the [Stan wiki](https://github.com/stan-dev/stan/wiki/Stan-Licensing).

Building tests with CMake
-------------------------

`runTests.py` uses CMake by default. A directory builds and runs one executable
containing all its tests; a source file builds and runs its individual executable:

```sh
./runTests.py -j6 ./test/unit/math/prim/fun
./runTests.py -j6 test/unit/math/prim/fun/sum_test.cpp
./runTests.py -j6 --build-dir .build -f sum test/unit/math/prim/fun
./runTests.py -j6 --cmake-arg=-DSTAN_THREADS=ON test/unit/math/rev/core
```

The first command selects `unit_math_prim_fun_all_tests`. `-f` filters the sources
compiled into the directory executable, and repeated filters must all match.
`--make-only` builds without running; `--test-only` runs the executable already
built with the same selection. `--changed` selects individual changed tests.
`--run-all` continues after failures and still returns a failure status.
`--jumbo-test` is accepted for compatibility; CMake directory builds already
combine the tests without generating amalgamated source files.
Build jobs are capped at six. XML results are written under
`<build-dir>/test-results`, and `build`, `.build`, and their suffixed variants
are ignored. CMake options go in `--cmake-arg`; `make/local` is used only by
`--backend make`. The existing Windows/ARM GitHub jobs explicitly retain Make.

The CMake test setup is adapted from `feature/cmake-tests` at
`5d304f35b23d355de9cc5958cbdc0df2abedb45f`. It uses the vendored dependency
versions in `make/libraries` and needs no downloads. Configure and build out
of source, with at most six compile jobs:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 6
ctest --test-dir build --parallel 6 --output-on-failure
```

Each test source has its own executable and CTest entry. Recursive build
targets such as `unit_math_rev_core_subtests` and `prob_math_bernoulli_subtests`
build a directory's tests. Select the matching CTest entries with `-R`, for
example `ctest --test-dir build -R '^test_unit_math_rev_core_' --output-on-failure`.
The optional `*_all_tests` targets link a directory's sources together; these
are excluded from the default build and CTest to avoid repeated execution.
Each test source compiles once into a shared object target, reused by its
individual executable and all containing directory executables. Switching
between a file, directory, and ancestor directory does not recompile that
test source. Tests use the repository root as their working directory.

Set `STAN_TEST_FILTER` to a regex matching repository-relative test source
paths to configure a smaller suite. For example:

```sh
cmake -S . -B build-core -DCMAKE_BUILD_TYPE=Release \
  -DSTAN_BUILD_PROB_TESTS=OFF -DSTAN_TEST_HEADERS=OFF \
  -DSTAN_TEST_FILTER='test/unit/math/(prim|fwd|rev|mix)/core/'
cmake --build build-core --parallel 6
ctest --test-dir build-core --parallel 6 --output-on-failure
```

`STAN_TEST_PCH` defaults to `ON`. Unit tests reuse a PCH for their autodiff
layer, or the AD helpers when the source already includes `test_ad.hpp`.
Grouped executables retain each source's PCH choice through shared objects.
Sources that change Stan or Eigen configuration macros skip PCH individually;
other sources in the same directory continue using their PCH. Umbrella and
multiple-translation-unit checks also compile without PCH.
Distribution tests reuse their fixture headers, without generated typedefs
in the PCH. Producers and consumers share compile flags and definitions;
GCC treats an invalid PCH as an error. Tests and vendored SUNDIALS share
`STAN_TEST_OPTIMIZATION` (default `1`, matching Make builds with `O=1`). Use
`-DSTAN_TEST_PCH=OFF` to check ordinary compilation. Standalone header
targets (`test-headers`, or a target
such as `stan_math_prim_fun_sum-test`) never use a PCH.

Distribution sources are generated during the build under `build/test/prob`,
using the existing `test/prob/generate_tests.cpp`. `STAN_PROB_TEST_BATCH_SIZE`
defaults to 100 cases per source. Each fixture's generation runs once before
its consumers and reruns when the fixture or generator changes.

`STAN_THREADS`, `STAN_OPENCL`, `STAN_MPI`, and `STAN_NO_RANGE_CHECKS` retain
their usual meanings. OpenCL and MPI tests are omitted when those features
are disabled. System dependencies can be selected with
`STAN_USE_SYSTEM_{BOOST,EIGEN,TBB,SUNDIALS}`; individual vendored roots can be
overridden with `STAN_{BOOST,EIGEN,TBB,SUNDIALS,OPENCL,BENCHMARK}_ROOT`.
The vendored TBB build currently supports GCC and Clang on Linux/macOS;
Windows builds require an installed TBB. MPI builds its Boost libraries from
the vendored sources, or uses installed libraries with `STAN_USE_SYSTEM_BOOST`. Set `BUILD_TESTING=OFF` to use the `stanmath`
interface target without configuring tests.

Jenkins uses the CMake/CTest pipeline with separate build directories per
configuration and no more than six compile jobs per build. ASan, GPU, MPI,
threading, and expression test configurations retain their coverage.
`STAN_TEST_ROW_VECTORS` and `STAN_PROB_TEST_ALL` enable the additional
distribution cases requested by Jenkins. Expression sources are prepared with
`python3 test/cmake_ci.py generate-expressions build-expressions`, then enabled
with `STAN_BUILD_EXPRESSION_TESTS=ON`. The runner also generates expression
sources when passed `test/expressions`, honoring `-e` and `--only-functions`.
Expression framework self-tests use the CMake runner. Var-matrix compatibility
compilation and linting retain Make.
