#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_ADJOINT_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_ADJOINT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/conditional_copy_and_promote.hpp>
#include <stan/math/rev/functor/reverse_pass_collect_adjoints.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <cstddef>
#include <ostream>
#include <tuple>
#include <utility>

namespace stan {
namespace math {
namespace internal {

/**
 * Build the reverse-mode result of a one-dimensional adaptive quadrature.
 *
 * The quadrature routine is supplied as a callable `integrator` so the three
 * `integrate_1d*` variants (plain, double-exponential, Gauss-Kronrod) can share
 * this body, differing only in which Boost routine and tolerances they bind.
 * `integrator` is invoked as `integrator(g)` where `g` is an integrand with
 * signature `g(x, xc)`; it returns the scalar integral of `g` over `[a, b]`.
 *
 * The integral value is computed with an all-`double` integrand. Adjoints are
 * then accumulated into the original `var` inputs:
 *  - endpoints contribute `-f(a)` and `f(b)` (skipped when infinite),
 *  - each `var` argument component contributes the integral of `d f / d arg`,
 *    obtained with nested reverse-mode autodiff.
 *
 * Gradients of `f` that come back NaN at a point where `f == 0` are set to
 * zero. Any other NaN propagates into the component integral and is reported as
 * a `domain_error` naming the (flattened) parameter index.
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Integrator callable `(integrand) -> double`
 * @tparam Args types of parameter pack arguments
 * @param function name of the calling function, used in error messages
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param integrator binds the quadrature routine and its tolerances
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename Integrator,
          typename... Args>
inline return_type_t<T_a, T_b, Args...> integrate_1d_adjoint(
    const char* function, const F& f, const T_a& a, const T_b& b,
    Integrator&& integrator, std::ostream* msgs, const Args&... args) {
  const double a_val = value_of(a);
  const double b_val = value_of(b);
  if (unlikely(a_val == b_val)) {
    if (is_inf(a_val)) {
      throw_domain_error(function, "Integration endpoints are both", a_val, "",
                         "");
    }
    return var(0.0);
  }
  auto args_val_tuple = std::make_tuple(value_of(args)...);
  auto eval_f = [&](const auto& x, const auto& xc) {
    return stan::math::apply(
        [&](auto&&... val_args) { return f(x, xc, msgs, val_args...); },
        args_val_tuple);
  };

  const double integral = integrator(eval_f);
  var ret(integral);
  if constexpr (is_var_v<T_a>) {
    double partial_a = 0.0;
    if (!is_inf(a_val)) {
      partial_a = -eval_f(a_val, 0.0);
    }
    reverse_pass_callback(
        [ret, a, partial_a]() { a.adj() += partial_a * ret.adj(); });
  }
  if constexpr (is_var_v<T_b>) {
    double partial_b = 0.0;
    if (!is_inf(b_val)) {
      partial_b = eval_f(b_val, 0.0);
    }
    reverse_pass_callback(
        [ret, b, partial_b]() { b.adj() += partial_b * ret.adj(); });
  }

  // Argument adjoints.
  if constexpr (is_any_var_scalar_v<Args...>) {
    auto args_adj = make_zeroed_arena(std::forward_as_tuple(args...));
    {
      nested_rev_autodiff argument_nest;
      auto args_copy = deep_copy_vargs<var>(std::forward_as_tuple(args...));
      auto args_copy_filter = filter_var_scalar_types(args_copy);
      auto integrate_grad = [&](auto&& target) -> double {
        return integrator([&](const auto& x, const auto& xc) {
          argument_nest.set_zero_all_adjoints();
          nested_rev_autodiff gradient_nest;
          var fx = stan::math::apply(
              [&](auto&&... local_args) {
                return f(x, xc, msgs, local_args...);
              },
              args_copy);
          fx.grad();
          double gradient = target.adj();
          if (is_nan(gradient) && fx.val() == 0) {
            gradient = 0.0;
          }
          return gradient;
        });
      };
      std::size_t param_index = 0;
      auto assign_grad = [&](auto&& adj, auto&& target) {
        adj = integrate_grad(target);
        if (unlikely(is_nan(adj))) {
          throw_domain_error("gradient_of_f", "The gradient of f", param_index,
                             "is nan for parameter ", "");
        }
        ++param_index;
      };
      iter_tuple_nested(
          [&](auto&& adj_leaf, auto&& var_leaf) {
            using leaf_t = std::decay_t<decltype(var_leaf)>;
            if constexpr (is_std_vector_v<leaf_t>) {
              for (size_t j = 0; j < var_leaf.size(); ++j) {
                assign_grad(adj_leaf[j], var_leaf[j]);
              }
            } else if constexpr (is_eigen_v<leaf_t>) {
              for (Eigen::Index j = 0; j < var_leaf.size(); ++j) {
                assign_grad(adj_leaf.coeffRef(j), var_leaf.coeff(j));
              }
            } else {  // scalar var leaf
              assign_grad(adj_leaf, var_leaf);
            }
          },
          args_adj, args_copy_filter);
    }
    auto args_filter = filter_var_scalar_types(std::forward_as_tuple(args...));
    reverse_pass_collect_adjoints(ret, args_filter, std::move(args_adj));
  }
  return ret;
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
