# Idiom Guide: I Need X, Use Y {#idioms}

This page maps common needs to the existing Stan Math symbol that already
covers them. Search it before writing a new helper, trait or loop. Paths are
relative to `stan/math/`. For the exact signatures of every symbol, generate the
API catalog with `./runClangd.py catalog` and grep `.agents/catalog/`.

### Types and metaprogramming (`prim/meta`)

| Need | Use | Header |
|---|---|---|
| Restrict an overload to certain types | `require_*_t<T>* = nullptr`. The `_vt` suffix tests the value type, `_st` the scalar type (see @ref require_meta_doc) | `prim/meta/is_*.hpp`, `require_generics.hpp` |
| Scalar type of a function's result | `return_type_t<T...>` | `prim/meta/return_type.hpp` |
| `double`-based type to hold partials | `partials_return_t<T...>` | `prim/meta/partials_return_type.hpp` |
| Innermost scalar of a container | `scalar_type_t<T>` | `prim/meta/scalar_type.hpp` |
| Element type of a container | `value_type_t<T>` | `prim/meta/value_type.hpp` |
| Arithmetic base type (strip all autodiff) | `base_type_t<T>` | `prim/meta/base_type.hpp` |
| Evaluated type of an Eigen expression | `plain_type_t<T>` | `prim/meta/plain_type.hpp` |
| Same container with a different scalar | `promote_scalar_t<S, T>` | `prim/meta/promote_scalar_type.hpp` |
| Is a type (or any of several) autodiff? | `is_autodiff_v<T>`, `is_any_autodiff_v<T...>` | `prim/meta/is_autodiff.hpp` |
| Is it `var_value<Matrix>`? | `is_var_matrix`, `require_var_matrix_t` | `prim/meta/is_var_matrix.hpp` |
| Return `var_value<Matrix>` or `Matrix<var>` to match the inputs | `conditional_var_value_t` | `rev/meta/conditional_var_value.hpp` |
| Skip constant terms when `propto` | `include_summand<propto, T...>` | `prim/meta/include_summand.hpp` |
| Return an expression that owns its temporaries | `make_holder` | `prim/meta/holder.hpp` |
| Branch on type at compile time | `if constexpr`. `forward_as` is legacy | |

### Functions (`prim/fun`)

| Need | Use | Header |
|---|---|---|
| Evaluate an expression used more than once | `to_ref`, `to_ref_if<cond>` | `prim/fun/to_ref.hpp` |
| Type that `to_ref` returns | `ref_type_t`, `ref_type_if_not_constant_t` | `prim/meta/ref_type.hpp` |
| Strip one level of autodiff | `value_of` | `prim/fun/value_of.hpp` |
| Strip all autodiff down to `double` | `value_of_rec` | `prim/fun/value_of_rec.hpp` |
| Index scalars and containers the same way | `scalar_seq_view`, `vector_seq_view` | `prim/fun/scalar_seq_view.hpp` |
| Values as an array or scalar (distributions) | `as_value_column_array_or_scalar`, `as_array_or_scalar`, `as_column_vector_or_scalar` | `prim/fun/` |
| Broadcast size, or empty check | `max_size`, `size_zero`, `math::size` | `prim/fun/` |
| Mathematical constants | `NEG_LOG_SQRT_TWO_PI`, `LOG_TWO`, `INFTY`, `NOT_A_NUMBER`, `EPSILON` | `prim/fun/constants.hpp` |
| Stable log-space arithmetic | `log_sum_exp`, `log1p_exp`, `log1m_exp`, `log_diff_exp` | `prim/fun/` |
| `x * log(y)` with `0 * log(0) = 0` | `lmultiply` | `prim/fun/lmultiply.hpp` |
| Force evaluation of an expression | `eval` | `prim/fun/eval.hpp` |

### Functors (`prim/functor`)

| Need | Use | Header |
|---|---|---|
| Apply a scalar function elementwise | `apply_scalar_unary<foo_fun, T>` with a `foo_fun::fun` struct (see `prim/fun/exp.hpp`) | `prim/functor/apply_scalar_unary.hpp` |
| Container-level math for arithmetic inputs | `apply_vector_unary<T>::apply`, `::reduce` | `prim/functor/apply_vector_unary.hpp` |
| Binary function with broadcasting | `apply_scalar_binary` | `prim/functor/apply_scalar_binary.hpp` |
| Hand-coded partial derivatives | `make_partials_propagator`, `partials<I>`, `partials_vec<I>`, `edge<I>` | `prim/functor/partials_propagator.hpp` |
| Loop over the elements of a tuple | `for_each` | `prim/functor/for_each.hpp` |
| Parallel sum over slices | `reduce_sum` | `prim/functor/reduce_sum.hpp` |
| Integrators and solvers | `integrate_1d`, `ode_rk45`, `map_rect`, `solve_newton`, `solve_powell` | `prim/functor/`, `rev/functor/` |

### Error checking (`prim/err`)

Checks take `(function, name, value)` and throw `std::domain_error` or
`std::invalid_argument` with a formatted message.

| Need | Use |
|---|---|
| Value checks | `check_positive`, `check_finite`, `check_not_nan`, `check_positive_finite`, `check_bounded` |
| Size and shape checks | `check_consistent_sizes`, `check_size_match`, `check_matching_dims`, `check_nonzero_size` |
| Matrix structure checks | `check_simplex`, `check_cholesky_factor`, `check_pos_definite` |
| A new elementwise check | `elementwise_check` (`prim/err/elementwise_check.hpp`) |
| Throw directly | `throw_domain_error` |

### Constraints (`prim/constraint`)

| Need | Use |
|---|---|
| Transform a value onto a constrained space, optionally with the Jacobian term | `lb_constrain`, `lub_constrain`, `simplex_constrain`, ... Overloads that take `lp` add the log Jacobian to it |
| The inverse transform | `lb_free`, `lub_free`, `simplex_free`, ... |

### Reverse mode (`rev`)

| Need | Use | Header |
|---|---|---|
| One output whose gradient is written as a closure | `make_callback_var`, `make_callback_vari` | `rev/core/callback_vari.hpp` |
| Custom adjoint code for several or matrix outputs | `reverse_pass_callback` | `rev/core/reverse_pass_callback.hpp` |
| Memory that survives until the reverse pass | `arena_t<T>`, `to_arena`, `arena_matrix`, `make_zeroed_arena` | `rev/meta/arena_type.hpp`, `rev/fun/to_arena.hpp`, `rev/core/` |
| Keep an object with a destructor for the reverse pass | `make_chainable_ptr` | `rev/core/chainable_object.hpp` |
| Gradients computed in advance | `precomputed_gradients` | `rev/core/precomputed_gradients.hpp` |
| Convert between `Matrix<var>` and `var_value<Matrix>` | `to_var_value`, `from_var_value` | `rev/fun/` |
| A nested gradient computation | `nested_rev_autodiff` | `rev/core/nested_rev_autodiff.hpp` |
| Gradient and Hessian drivers | `gradient`, `hessian` | `rev/functor/gradient.hpp`, `mix/functor/hessian.hpp` |

### Forward mode (`fwd`)

- An `fvar` specialization returns `fvar<T>(f(x.val_), x.d_ * df(x.val_))`.
  See `fwd/fun/exp.hpp`.

### Tests (`test/unit/math`)

| Need | Use |
|---|---|
| Test values and gradients in all modes | `stan::test::expect_ad` (`test_ad.hpp`) |
| Elementwise functions | `expect_ad_vectorized`, `expect_ad_vectorized_binary` |
| `var_value<Matrix>` support | `expect_ad_matvar` (`test_ad_matvar.hpp`) |
| Every mode throws | `expect_all_throw` |
| Custom tolerances | `ad_tolerances` |
| Relative near-equality | `expect_near_rel` |
| Rev test fixture | `TEST_F(AgradRev, ...)` (`rev/util.hpp`) |
