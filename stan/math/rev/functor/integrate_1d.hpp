#ifndef STAN_MATH_REV_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_REV_FUNCTOR_integrate_1d_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/core/precomputed_gradients.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/integrate_1d.hpp>
#include <cmath>
#include <functional>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the integral of f from a to b to the given relative tolerance
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
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_var<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &... args) {
  static constexpr const char *function = "integrate_1d";
  check_less_or_equal(function, "lower limit", a, b);

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

    double integral = integrate(
        [&](const auto &x, const auto &xc) {
          return apply(
              [&](auto &&... val_args) { return f(x, xc, msgs, val_args...); },
              args_val_tuple);
        },
        a_val, b_val, relative_tolerance);

    constexpr size_t num_vars_ab = is_var<T_a>::value + is_var<T_b>::value;
    size_t num_vars_args = count_vars(args...);
    vari **varis = ChainableStack::instance_->memalloc_.alloc_array<vari *>(
        num_vars_ab + num_vars_args);
    double *partials = ChainableStack::instance_->memalloc_.alloc_array<double>(
        num_vars_ab + num_vars_args);
    // We move this pointer up based on whether we a or b is a var type.
    double *partials_ptr = partials;

    save_varis(varis, a, b, args...);

    for (size_t i = 0; i < num_vars_ab + num_vars_args; ++i) {
      partials[i] = 0.0;
    }

    if (is_var<T_a>::value && !is_inf(a)) {
      *partials_ptr = apply(
          [&f, a_val, msgs](auto &&... val_args) {
            return -f(a_val, 0.0, msgs, val_args...);
          },
          args_val_tuple);
      partials_ptr++;
    }

    if (!is_inf(b) && is_var<T_b>::value) {
      *partials_ptr = apply(
          [&f, b_val, msgs](auto &&... val_args) {
            return f(b_val, 0.0, msgs, val_args...);
          },
          args_val_tuple);
      partials_ptr++;
    }

    {
      nested_rev_autodiff argument_nest;
      // The arguments copy is used multiple times in the following nests, so
      // do it once in a separate nest for efficiency
      auto args_tuple_local_copy = std::make_tuple(deep_copy_vars(args)...);

      // Save the varis so it's easy to efficiently access the nth adjoint
      std::vector<vari *> local_varis(num_vars_args);
      apply(
          [&](const auto &... args) {
            save_varis(local_varis.data(), args...);
          },
          args_tuple_local_copy);

      for (size_t n = 0; n < num_vars_args; ++n) {
        // This computes the integral of the gradient of f with respect to the
        // nth parameter in args using a nested nested reverse mode autodiff
        *partials_ptr = integrate(
            [&](const auto &x, const auto &xc) {
              argument_nest.set_zero_all_adjoints();

              nested_rev_autodiff gradient_nest;
              var fx = apply(
                  [&f, &x, &xc, msgs](auto &&... local_args) {
                    return f(x, xc, msgs, local_args...);
                  },
                  args_tuple_local_copy);
              fx.grad();

              double gradient = local_varis[n]->adj();

              // Gradients that evaluate to NaN are set to zero if the function
              // itself evaluates to zero. If the function is not zero and the
              // gradient evaluates to NaN, a std::domain_error is thrown
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
            a_val, b_val, relative_tolerance);
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
 * Compute the integral of the single variable function f from a to b to within
 * a specified relative tolerance. a and b can be finite or infinite.
 *
 * f should be compatible with reverse mode autodiff and have the signature:
 *   var f(double x, double xc, const std::vector<var>& theta,
 *     const std::vector<double>& x_r, const std::vector<int> &x_i,
 * std::ostream* msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream.
 *
 * Integrals that cross zero are broken into two, and the separate integrals are
 * each integrated to the given relative tolerance.
 *
 * For integrals with finite limits, the xc argument is the distance to the
 * nearest boundary. So for a > 0, b > 0, it will be a - x for x closer to a,
 * and b - x for x closer to b. xc is computed in a way that avoids the
 * precision loss of computing a - x or b - x manually. For integrals that cross
 * zero, xc can take values a - x, -x, or b - x depending on which integration
 * limit it is nearest.
 *
 * If either limit is infinite, xc is set to NaN
 *
 * The integration algorithm terminates when
 *   \f[
 *     \frac{{|I_{n + 1} - I_n|}}{{|I|_{n + 1}}} < \text{relative tolerance}
 *   \f]
 * where \f$I_{n}\f$ is the nth estimate of the integral and \f$|I|_{n}\f$ is
 * the nth estimate of the norm of the integral.
 *
 * Integrals that cross zero are
 * split into two. In this case, each integral is separately integrated to the
 * given relative_tolerance.
 *
 * Gradients of f that evaluate to NaN when the function evaluates to zero are
 * set to zero themselves. This is due to the autodiff easily overflowing to NaN
 * when evaluating gradients near the maximum and minimum floating point values
 * (where the function should be zero anyway for the integral to exist)
 *
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam T_theta type of parameters
 * @tparam T Type of f
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename T_theta,
          typename = require_any_var_t<T_a, T_b, T_theta>>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance = std::sqrt(EPSILON)) {
  return integrate_1d_impl(integrate_1d_adapter<F>(f), a, b, relative_tolerance,
                           msgs, theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan

#endif
