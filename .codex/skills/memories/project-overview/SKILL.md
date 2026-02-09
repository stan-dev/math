# Project Overview

- Purpose: Stan Math Library is a C++ reverse-mode automatic differentiation library used by Stan (see `README.md`).
- Core layout: primary headers and implementation live under `stan/math/` with domain splits in `stan/math/prim`, `stan/math/rev`, `stan/math/fwd`, `stan/math/mix`, `stan/math/opencl`. Umbrella headers include `stan/math.hpp`, `stan/math/prim.hpp`, `stan/math/rev.hpp`, `stan/math/fwd.hpp`, `stan/math/mix.hpp`.
- Tests: unit tests are under `test/unit/` (notably `test/unit/math/...`); expression tests under `test/expressions/`; probability tests and generators under `test/prob/`.
- Benchmarks: microbenchmarks live in `benchmarks/`.
- Third-party deps: vendored libraries are under `lib/` (TBB, Sundials, etc. referenced in build rules).

Golden paths
- Add a function: implement in the appropriate domain folder (e.g., `stan/math/prim/...`), update umbrella headers if needed, and add a unit test under `test/unit/math/...`.
- Fix a bug: find the failing test under `test/unit/...`, adjust code in `stan/math/...`, and run the focused test via `./runTests.py` or `make` (see `makefile`).
- Performance work: add or update a benchmark in `benchmarks/` and compare results locally.

Risk/fragile areas
- OpenCL code paths (`stan/math/opencl/`) and SUNDIALS-dependent tests may require extra dependencies (see `make/tests`).
