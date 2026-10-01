#ifndef STAN_MATH_PRIM_FUN_PHI_HPP
#define STAN_MATH_PRIM_FUN_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * The unit normal cumulative distribution function.
 *
 * The return value for a specified input is the probability that
 * a random unit normal variate is less than or equal to the
 * specified value, defined by
 *
 * \f$\Phi(x) = \int_{-\infty}^x \mbox{\sf Norm}(x|0, 1) \ dx\f$
 *
 * This function can be used to implement the inverse link function
 * for probit regression.
 *
 * @param x Argument.
 * @return Probability random sample is less than or equal to argument.
 */
inline double Phi(double x) {
  check_not_nan("Phi", "x", x);
  return exp(internal::std_normal_lcdf_value_grad<false>(x).first);
}

/**
 * Structure to wrap Phi() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Unit normal CDF of x.
 */
struct Phi_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return Phi(std::forward<T>(x));
  }
};

/**
 * Vectorized version of Phi().
 *
 * @tparam T type of container
 * @param x container
 * @return Unit normal CDF of each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto Phi(T&& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
