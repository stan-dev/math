# Tests

Guide: `doxygen/contributor_help_pages/autodiff_test_guide.md`.

- **Layout.** `test/unit/math/{prim,rev,fwd,mix,opencl}/<area>/foo_test.cpp`
  mirrors `stan/math/`. Every file must end in `_test.cpp`, and gtest names
  must be unique across the repo (`runChecks.py`).
- **Mix tests (primary).** Mix tests cover gradients for every autodiff
  type.
  - Include `test/unit/math/test_ad.hpp`.
  - Write a lambda `f` over `auto` arguments and call
    `stan::test::expect_ad(f, x...)` with one to three arguments.
  - Elementwise functions: `expect_ad_vectorized` /
    `expect_ad_vectorized_binary`.
  - `var_value<Matrix>`: `expect_ad_matvar`
    (`test/unit/math/test_ad_matvar.hpp`).
  - Tolerances: pass `stan::test::ad_tolerances` as the first argument.
  - Invalid inputs: `expect_ad` already checks that every mode throws
    consistently. Use `expect_all_throw` for throw-only checks.
  - Example: `test/unit/math/mix/fun/exp_test.cpp`.
- **Rev tests.** Use `TEST_F(AgradRev, Name)` from
  `test/unit/math/rev/util.hpp`; never plain `TEST(`.
- **Prim tests.** Include `<stan/math/prim.hpp>` and use doubles only.
- **Distributions.** `test/prob/` is generated (`make generate-tests`). Add
  a `test/prob/<dist>/<dist>_test.hpp` there rather than hand-writing
  gradient tests.
- **Running.** `./runTests.py <file>`; `-f <substring>` filters by path.
  Save long output to a file and grep it; don't pipe it through `tail`.
