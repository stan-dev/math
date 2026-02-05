# Style and Conventions

Formatting
- C++ formatting follows `.clang-format` (Google-based) with `ColumnLimit: 80`, `IndentWidth: 2`, `PointerAlignment: Left`, and `SortIncludes: false`.

Code organization
- Library code lives in `stan/math/` with domains in `prim`, `rev`, `fwd`, `mix`, and `opencl`.
- Public umbrella headers exist at `stan/math.hpp` and per-domain headers like `stan/math/prim.hpp`.

Naming and file patterns
- Headers use `.hpp`.
- Unit tests are `*_test.cpp` under `test/unit/...`.
- Header compile tests use the `-test` suffix (e.g., `stan/math/constants.hpp-test`).

Quality checks
- Use `make cpplint` and `make clang-tidy` for static checks (targets defined in `make/cpplint` and `make/clang-tidy`).
