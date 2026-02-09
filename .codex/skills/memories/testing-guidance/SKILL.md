# Testing Guidance

- Primary runner: `./runTests.py` (see help with `-h`).
- Unit tests live under `test/unit/`; run a subtree with `./runTests.py test/unit` or a specific test file path.
- Build a single test binary with `make test/unit/path/to/foo_test` and run the produced binary.
- Header compile tests: `make test-headers` or `make stan/math/constants.hpp-test`.
- Dependency checks: `./runChecks.py` or `make test-math-dependencies`.
