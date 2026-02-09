# Suggested Commands

Setup/overview
- `make help`: print available make targets and testing guidance (see `makefile`).

Testing
- `./runTests.py -h`: show test runner help and options.
- `./runTests.py test/unit`: run unit tests under `test/unit/`.
- `./runTests.py -f <pattern>`: run tests matching a pattern (e.g., `-f opencl`).
- `./runTests.py --changed`: run tests for files changed in git.
- `make test/unit/math/prim/fun/abs_test`: build a single unit test binary.
- `make test-headers`: compile-check all headers (see `make/tests`).
- `make stan/math/constants.hpp-test`: compile-check a single header.

Quality
- `make cpplint`: run cpplint over source files.
- `make clang-tidy`: run clang-tidy across tests (optional `files=...`).
- `make clang-tidy-fix`: clang-tidy with autofix output to `.clang-fixes.yml`.

Docs
- `make doxygen`: generate API docs in `doc/api/`.

Cleanup
- `make clean`: remove generated test files/binaries.
- `make clean-deps`: remove dependency files if tests stop building.
- `make clean-all`: deep clean including docs/libraries.

Utilities
- `./runChecks.py`: check math dependency rules (targeted by `make test-math-dependencies`).
