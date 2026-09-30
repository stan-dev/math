# Stan Math conventions for new or extended code

## Enforced by CI

`make test-math-dependencies` runs `./runChecks.py`. It greps the source,
skipping comments, and fails on each of the following.

**Files under `stan/math/prim`**
- Including `<stan/math/rev/`, `<stan/math/fwd/` or `<stan/math/mix/`.
- Naming `stan::math::var` or `stan::math::fvar`.

**Files under `stan/math/rev`**
- Including `<stan/math/fwd/` or `<stan/math/mix/`.
- Naming `stan::math::fvar`.

**Anywhere in `stan/math`**
- Boost type traits and enable_if:
  - `boost::is_unsigned`, `is_arithmetic`, `is_convertible`, `is_same`
  - `boost::enable_if`, `enable_if_c`, `disable_if`
  - their `<boost/...>` headers
  
  Use `<type_traits>` instead. In practice, prefer Stan's `require_*` aliases
  over a raw `std::enable_if`.
- `std::lgamma`, which is not reentrant. Use Stan's `lgamma`.

**Tests**
- Every `.cpp` under `test/unit` must end in `_test.cpp`. A `*_test.hpp` there
  is an error.
- Tests in `test/unit/math/rev` must use `TEST_F(AgradRev, ...)` (from
  `test/unit/math/rev/util.hpp`) or another fixture form, not a raw `TEST(`.
- `TEST`/`TEST_F` names must be unique across `test/unit`.

## Layout

- A function `foo` lives in `stan/math/prim/fun/foo.hpp`, which handles
  arithmetic types and containers.
- Autodiff specializations go in `stan/math/rev/fun/foo.hpp` and
  `stan/math/fwd/fun/foo.hpp`.
- Include each new header from the matching aggregate header:
  `stan/math/prim/fun.hpp`, `rev/fun.hpp`, `fwd/fun.hpp`, and the same for
  `meta.hpp`, `err.hpp`, and so on.
- Every header must compile on its own:
  `make stan/math/prim/fun/foo.hpp-test`.
- Tests: `test/unit/math/mix/fun/foo_test.cpp` using `stan::test::expect_ad`
  (see `doxygen/contributor_help_pages/autodiff_test_guide.md`).

## Public entry points

- Validate arguments with the `check_*` functions from `prim/err`, whose
  signature is `(function, name, y)`. Declare
  `static constexpr const char* function = "foo";` once.
- Constrain overloads with the `require_*` aliases (template parameter
  `require_..._t<T>* = nullptr`), not hand-written SFINAE.
- Take Eigen arguments as generic templates so expressions are accepted, not
  as `const Eigen::MatrixXd&`.
- If an Eigen argument is read more than once, evaluate it once with
  `to_ref(x)` (type `ref_type_t<T>`).
- Returning an expression that refers to a local temporary needs
  `make_holder`.

## Reverse mode (rev)

- Prefer `make_callback_var(value, [captures](auto& vi) mutable {...})`
  (`rev/core/callback_vari.hpp`) or `reverse_pass_callback`
  (`rev/core/reverse_pass_callback.hpp`) over a new `vari` subclass.
- Everything a reverse-pass lambda captures must be a `var`, `var_value<T>`,
  or arena-backed storage (`arena_t<T>`, `to_arena(x)`).
  - Capturing an `Eigen::Matrix`, a `std::vector`, or a reference to a local
    object dangles or leaks, because arena memory never runs destructors.
- Write matrix overloads for both representations:
  - `var_value<Eigen::Matrix>`, via `require_var_matrix_t`
  - `Eigen::Matrix<var>`, via `require_eigen_vt<is_var, T>`
  
  See `doxygen/contributor_help_pages/reverse_mode_types.md`.
- Never hold an Eigen expression in an `auto` variable in rev code. Evaluate it
  or use a `to_ref`.
- Don't return non-arena temporaries from code that is captured for the
  reverse pass.

## Compile-time branching

- Use `if constexpr (is_autodiff_v<T>)` (or `is_any_autodiff_v<Ts...>`) to skip
  gradient code for constants.
- `forward_as` is legacy and unused (only stale includes remain);
  `is_constant_all` survives in older code only. Don't use either.

## Distributions (prim/prob)

- Follow `stan/math/prim/prob/normal_lpdf.hpp`:
  - forwarding-reference arguments
  - `ref_type_if_not_constant_t`
  - `check_consistent_sizes`, then the value checks
  - early returns on `size_zero` and `!include_summand<propto, ...>`
  - `make_partials_propagator`, `partials<I>(ops)` (0-based), `ops.build(logp)`
- Don't use `operands_and_partials`/`edgeN_`, and don't add new `*_log`
  aliases.

## Formatting

- The repo uses clang-format 10 through `hooks/pre-commit`. A newer
  clang-format rewraps unrelated code, so format only the lines you changed.
- `make cpplint` must pass.
