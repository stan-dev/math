#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/integrate_1d_adjoint.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
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
    const Args &... args) {
  static constexpr const char *function = "integrate_1d_double_exponential";
  check_less_or_equal(function, "lower limit", a, b);
  check_nonnegative(function, "max_refinements", max_refinements);
  check_nonnegative(function, "absolute_tolerance", absolute_tolerance);
  return internal::integrate_1d_adjoint(
      function, f, a, b,
      [&](auto &&integrand) {
        return integrate_de(std::forward<decltype(integrand)>(integrand),
                            value_of(a), value_of(b), relative_tolerance,
                            absolute_tolerance, max_refinements);
      },
      msgs, args...);
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
    const Args &... args) {
  return integrate_1d_double_exponential_tol(
      f, a, b, std::sqrt(EPSILON), 0.0,
      INTEGRATE_1D_DOUBLE_EXPONENTIAL_MAX_REFINEMENTS, msgs, args...);
}

}  // namespace math
}  // namespace stan

#endif
