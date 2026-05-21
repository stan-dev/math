#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/integrate_1d_double_exponential.hpp>
#include <cmath>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the integral of f from a to b using adaptive double-exponential
 * quadrature.
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Args types of parameter pack arguments
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the Boost
 *   quadrature class constructor
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_double_exponential_tol(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    double absolute_tolerance, int max_refinements, std::ostream *msgs,
    const Args &...args) {
  static constexpr const char *function = "integrate_1d_double_exponential";
  check_less_or_equal(function, "lower limit", a, b);
  check_nonnegative(function, "max_refinements", max_refinements);
  check_nonnegative(function, "absolute_tolerance", absolute_tolerance);

  double a_val = value_of(a);
  double b_val = value_of(b);

  if (unlikely(a_val == b_val)) {
    if (is_inf(a_val)) {
      throw_domain_error(function, "Integration endpoints are both", a_val, "",
                         "");
    }
    return var(0.0);
  } else {
    auto args_val_tuple = std::make_tuple(value_of(args)...);

    double integral = integrate_de(
        [&](const auto &x, const auto &xc) {
          return math::apply(
              [&](auto &&...val_args) { return f(x, xc, msgs, val_args...); },
              args_val_tuple);
        },
        a_val, b_val, relative_tolerance, absolute_tolerance, max_refinements);

    constexpr size_t num_vars_ab = is_var<T_a>::value + is_var<T_b>::value;
    size_t num_vars_args = count_vars(args...);
    vari **varis = ChainableStack::instance_->memalloc_.alloc_array<vari *>(
        num_vars_ab + num_vars_args);
    double *partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        num_vars_ab + num_vars_args);
    double *partials_ptr = partials;

    save_varis(varis, a, b, args...);

    for (size_t i = 0; i < num_vars_ab + num_vars_args; ++i) {
      partials[i] = 0.0;
    }

    if constexpr (is_var<T_a>::value) {
      if (!is_inf(a)) {
        *partials_ptr = math::apply(
            [&f, a_val, msgs](auto &&...val_args) {
              return -f(a_val, 0.0, msgs, val_args...);
            },
            args_val_tuple);
        partials_ptr++;
      }
    }

    if constexpr (is_var<T_b>::value) {
      if (!is_inf(b)) {
        *partials_ptr = math::apply(
            [&f, b_val, msgs](auto &&...val_args) {
              return f(b_val, 0.0, msgs, val_args...);
            },
            args_val_tuple);
        partials_ptr++;
      }
    }

    {
      nested_rev_autodiff argument_nest;
      auto args_tuple_local_copy = std::make_tuple(deep_copy_vars(args)...);
      std::vector<vari *> local_varis(num_vars_args);
      math::apply(
          [&](const auto &...args) { save_varis(local_varis.data(), args...); },
          args_tuple_local_copy);

      for (size_t n = 0; n < num_vars_args; ++n) {
        *partials_ptr = integrate_de(
            [&](const auto &x, const auto &xc) {
              argument_nest.set_zero_all_adjoints();

              nested_rev_autodiff gradient_nest;
              var fx = math::apply(
                  [&f, &x, &xc, msgs](auto &&...local_args) {
                    return f(x, xc, msgs, local_args...);
                  },
                  args_tuple_local_copy);
              fx.grad();

              double gradient = local_varis[n]->adj();

              if (is_nan(gradient)) {
                if (fx.val() == 0) {
                  gradient = 0;
                } else {
                  throw_domain_error("gradient_of_f", "The gradient of f", n,
                                     "is nan for parameter ", "");
                }
              }
              return gradient;
            },
            a_val, b_val, relative_tolerance, absolute_tolerance,
            max_refinements);
        partials_ptr++;
      }
    }

    return make_callback_var(
        integral,
        [total_vars = num_vars_ab + num_vars_args, varis, partials](auto &vi) {
          for (size_t i = 0; i < total_vars; ++i) {
            varis[i]->adj_ += partials[i] * vi.adj();
          }
        });
  }
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive double-exponential quadrature. a and b can be finite or
 * infinite.
 *
 * f should be compatible with reverse mode autodiff and have the signature:
 *   var f(double x, double xc, const std::vector<var>& theta,
 *     const std::vector<double>& x_r, const std::vector<int> &x_i,
 *     std::ostream* msgs)
 *
 * Gradients of f that evaluate to NaN when the function evaluates to zero
 * are set to zero themselves.
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Args types of parameter pack arguments
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the Boost
 *   quadrature class constructor
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_double_exponential(
    const F &f, const T_a &a, const T_b &b, std::ostream *msgs,
    const Args &...args) {
  return integrate_1d_double_exponential_tol(
      f, a, b, std::sqrt(EPSILON), 0.0,
      INTEGRATE_1D_DOUBLE_EXPONENTIAL_MAX_REFINEMENTS, msgs, args...);
}

}  // namespace math
}  // namespace stan

#endif
