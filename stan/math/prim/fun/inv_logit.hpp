#ifndef STAN_MATH_PRIM_FUN_INV_LOGIT_HPP
#define STAN_MATH_PRIM_FUN_INV_LOGIT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the inverse logit function applied to the argument.
 *
 * The inverse logit function is defined by
 *
 * \f$\mbox{logit}^{-1}(x) = \frac{1}{1 + \exp(-x)}\f$.
 *
 * This function can be used to implement the inverse link function
 * for logistic regression.
 *
 * The inverse to this function is <code>logit</code>.
 *
 \f[
 \mbox{inv\_logit}(y) =
 \begin{cases}
 \mbox{logit}^{-1}(y) & \mbox{if } -\infty\leq y \leq \infty \\[6pt]
 \textrm{NaN} & \mbox{if } y = \textrm{NaN}
 \end{cases}
 \f]

 \f[
 \frac{\partial\, \mbox{inv\_logit}(y)}{\partial y} =
 \begin{cases}
 \frac{\partial\, \mbox{logit}^{-1}(y)}{\partial y} & \mbox{if } -\infty\leq
 y\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } y = \textrm{NaN} \end{cases} \f]

 \f[
 \mbox{logit}^{-1}(y) = \frac{1}{1 + \exp(-y)}
 \f]

 \f[
 \frac{\partial \, \mbox{logit}^{-1}(y)}{\partial y} =
 \frac{\exp(y)}{(\exp(y)+1)^2} \f]
 *
 * @param a Argument.
 * @return Inverse logit of argument.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline double inv_logit(T&& a) {
  if (a < 0) {
    double exp_a = std::exp(a);
    if (a < LOG_EPSILON) {
      return exp_a;
    }
    return exp_a / (1.0 + exp_a);
  }
  return inv(1.0 + std::exp(-a));
}

/**
 * Structure to wrap inv_logit() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Inverse logit of x.
 */
struct inv_logit_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return inv_logit(std::forward<T>(x));
  }
};

/**
 * Vectorized version of inv_logit() for containers containing ad types.
 *
 * @tparam T type of std::vector
 * @param x std::vector
 * @return Inverse logit applied to each value in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr,
          require_not_rev_matrix_t<Container>* = nullptr>
inline auto inv_logit(Container&& x) {
  return apply_scalar_unary<inv_logit_fun, Container>::apply(
      std::forward<Container>(x));
}

/**
 * Vectorized version of inv_logit() for containers with arithmetic scalar
 * types.
 *
 * @tparam T A type of either `std::vector` or a type that directly inherits
 * from `Eigen::DenseBase`. The inner scalar type must not have an auto diff
 * scalar type.
 * @param x Eigen expression
 * @return Inverse logit applied to each value in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              Container>* = nullptr>
inline auto inv_logit(Container&& x) {
  return apply_vector_unary<Container>::apply(
      std::forward<Container>(x),
      [](auto&& v) { return v.array().logistic(); });
}
}  // namespace math
}  // namespace stan

#endif
