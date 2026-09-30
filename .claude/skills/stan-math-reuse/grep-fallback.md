# Grep fallback (no catalog)

Use these when `.agents/catalog/` is missing. Generate the catalog with
`./runClangd.py catalog`.

All recipes are scoped to `stan/math` and run from the repo root. Replace `NAME`
with the symbol or keyword. Add `| head` when output may be long.

## Where is a symbol defined?

File names usually match function names (`stan/math/prim/fun/exp.hpp` defines
`exp`). Some helpers live elsewhere:

| Symbol | Header |
|---|---|
| `make_holder` | `prim/meta/holder.hpp` |
| `make_callback_var` | `rev/core/callback_vari.hpp` |
| `arena_t` | `rev/meta/arena_type.hpp` |
| `to_arena` | `rev/fun/to_arena.hpp` |

Type aliases, traits and classes:

```sh
git grep -n -E "(using NAME\b\s*=|struct NAME\b|class NAME\b)" -- stan/math
```

Functions, every overload across prim/rev/fwd:

```sh
git grep -n -E "^\s*(inline\s+)?[^=;]*\bNAME\(" -- 'stan/math/*/fun/*.hpp' | grep inline
```

## Search by behavior (doc comments)

The first sentence of each `/** */` block is its brief.

```sh
git grep -n -i -E "^\s*\*\s.*KEYWORD" -- stan/math/prim/fun | head
ls stan/math/prim/fun | grep -i KEYWORD
```

For example, `ls stan/math/prim/fun | grep -i -E "log.*exp|exp.*log"` finds
`log_sum_exp`, `log1p_exp`, `log1m_exp` and `log_diff_exp`.

## Type traits and SFINAE (prim/meta)

- Before writing `enable_if` or a new trait, list the existing `require_*`
  aliases for a trait.
- The pattern is `require_[all_|any_|not_|all_not_|any_not_]<trait>_[t|vt|st]`:
  - `_vt` checks the `value_type`
  - `_st` checks the `scalar_type`

```sh
git grep -n -E "^\s*using (require_\w*TRAIT\w*) =" -- stan/math/prim/meta stan/math/rev/meta
git grep -n -E "^struct is_\w*TRAIT" -- stan/math/prim/meta stan/math/rev/meta
```

Common type helpers, all in `prim/meta/<name>.hpp`:
- `value_type_t`
- `scalar_type_t`
- `return_type_t`
- `partials_return_t`
- `plain_type_t`
- `ref_type_t`
- `promote_scalar_t` (`promote_scalar_type.hpp`)
- `base_type_t`

## Error checks (prim/err)

`check_*` functions take `(function, name, y)`. List them all with:

```sh
ls stan/math/prim/err
```

For a new elementwise check, see `elementwise_check` in
`prim/err/elementwise_check.hpp`.

## Vectorization and functors (prim/functor)

```sh
ls stan/math/prim/functor
git grep -n -E "struct \w+_fun\b" -- stan/math/prim/fun | head
```

- `apply_scalar_unary<foo_fun, T>` uses the `foo_fun::fun` pattern; see
  `prim/fun/exp.hpp`.
- `apply_vector_unary<T>::apply` / `::reduce` handles whole-container math.
- `apply_scalar_binary` handles two-argument broadcasting.

## Reverse mode (rev/core)

```sh
git grep -l "reverse_pass_callback\|make_callback_var" -- stan/math/rev/fun | head
```

Read one of those files as a template before writing a new `vari`.

## Existing callers of a candidate

```sh
git grep -l -w NAME -- stan/math | head
```
