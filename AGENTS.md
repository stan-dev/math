# Stan Math: agent guide

Header-only C++17 automatic differentiation library. Most code lives in
`stan/math/{prim,rev,fwd,mix,opencl}`; tests live in `test/unit/math/`.

## Reuse before you write
- The library already has 1900+ headers. Most helpers and traits you might write already exist.
- Before adding any function, overload, trait, `require_*` alias, struct or container loop under `stan/`, find what already exists:
  1. Use an LSP if available to query for existing symbols.
  1. Check the idiom guide, `doxygen/contributor_help_pages/idioms.md` ("I need X, use Y").
  2. Grep the API catalog in `.agents/catalog/` (one line per symbol: `name | header | signature | brief`). Generate the catalog with `./runClangd.py catalog`; it is gitignored and never committed.
  3. Grep `stan/math/` for real call sites of the candidates.
- File names usually equal the function name (`prim/fun/foo.hpp` defines `foo`). Exceptions include `make_holder` (`prim/meta/holder.hpp`), `make_callback_var` (`rev/core/callback_vari.hpp`) and `arena_t` (`rev/meta/arena_type.hpp`). Type traits that end in `_t` will usually be in a file missing the `_t` such as `scalar_type.hpp` which has `scalar_type_t`.

## Layers
- `prim` works for any scalar type, `rev` specializes for `var`
  (reverse mode), `fwd` for `fvar` (forward mode), and `mix` for nested types. `opencl` is GPU code behind `STAN_OPENCL`.
- `prim` must not include `rev`/`fwd`/`mix` or name `var`/`fvar`, and `rev` must not include `fwd`/`mix`. `make test-math-dependencies` (`runChecks.py`) enforces this.
- One function `foo` has `prim/fun/foo.hpp`, plus optional `rev/fun/foo.hpp` and `fwd/fun/foo.hpp` for specialized derivatives. Add each new header to the matching aggregate header (`prim/fun.hpp`, `rev/fun.hpp`, `fwd/fun.hpp`, `prim/prob.hpp`).
- Scalar functions are vectorized with `apply_scalar_unary` and `apply_vector_unary`; follow the pattern in `prim/fun/exp.hpp`.

## Commands
- Run a test: `./runTests.py test/unit/math/mix/fun/foo_test`. Useful flags: `-j N`, `-f <substring>`, `--changed`. Keep `-j` at less than half the available logical cores.
- Check a header includes what it uses: `make stan/math/prim/fun/foo.hpp-test` (all headers: `make test-headers`).
- Lint and layering: `make cpplint`, `make test-math-dependencies`.
- Local build options (for example `STAN_OPENCL=true`) go in `make/local`.
- Format with clang-format through `hooks/pre-commit` (install with `hooks/install_hooks.sh`). Newer clang-format versions rewrap unrelated code, so never reformat whole files.
- PR checklist: `.github/PULL_REQUEST_TEMPLATE.md`.

## Human docs (read the relevant one before starting)
All are in `doxygen/contributor_help_pages/`:
- `getting_started.md`: the full recipe for adding a function.
- `common_pitfalls.md`: Eigen expressions, `auto`, holders and the arena.
- `require_meta.md`: the `require_*` SFINAE traits.
- `reverse_mode_types.md`: `var_value<Matrix>` vs `Matrix<var>`.
- `adding_new_distributions.md`, `distribution_tests.md`: `_lpdf`/`_lpmf`.
- `autodiff_test_guide.md`: `expect_ad` testing.
- `add_new_opencl_kernel.md`: GPU code.

Some directories have their own `AGENTS.md` (`stan/math/prim/prob`,
`stan/math/rev`, `stan/math/opencl`, `test`).
