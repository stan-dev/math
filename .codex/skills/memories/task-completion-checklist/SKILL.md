# Task Completion Checklist

- [ ] Run focused tests for changed areas (e.g., `./runTests.py test/unit/...` or single test binary via `make`).
- [ ] Add/update tests when behavior changes; ensure new data files are referenced correctly.
- [ ] Run relevant linters/static checks if touching core code (`make cpplint`, `make clang-tidy`).
- [ ] Update docs/changelogs if user-facing behavior changes (rare for math internals).
- [ ] Self-review for dependency rules (prim/fwd/rev/mix layering) and header include correctness.
