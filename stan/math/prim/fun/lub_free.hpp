#ifndef STAN_MATH_PRIM_FUN_LUB_FREE_HPP
#define STAN_MATH_PRIM_FUN_LUB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/logit.hpp>
#include <stan/math/prim/fun/lb_free.hpp>
#include <stan/math/prim/fun/ub_free.hpp>

namespace stan {
namespace math {

/**
 * Return the unconstrained scalar that transforms to the
 * specified lower- and upper-bounded scalar given the specified
 * bounds.
 *
 * <p>The transform in <code>lub_constrain(T, double, double)</code>,
 * is reversed by a transformed and scaled logit,
 *
 * <p>\f$f^{-1}(y) = \mbox{logit}(\frac{y - L}{U - L})\f$
 *
 * where \f$U\f$ and \f$L\f$ are the lower and upper bounds.
 *
 * <p>If the lower bound is negative infinity and upper bound finite,
 * this function reduces to <code>ub_free(y, ub)</code>.  If
 * the upper bound is positive infinity and the lower bound
 * finite, this function reduces to
 * <code>lb_free(x, lb)</code>.  If the upper bound is
 * positive infinity and the lower bound negative infinity,
 * this function reduces to <code>identity_free(y)</code>.
 *
 * @tparam T type of scalar
 * @tparam L type of lower bound
 * @tparam U type of upper bound
 * @param y constrained value
 * @param lb lower bound
 * @param ub upper bound
 * @return the free scalar that transforms to the input scalar
 *   given the bounds
 * @throw std::invalid_argument if the lower bound is greater than
 *   the upper bound, y is less than the lower bound, or y is
 *   greater than the upper bound
 */
template <typename T, typename L, typename U,
          require_stan_scalar_t<T>* = nullptr>
inline auto lub_free(T&& y, const L& lb, const U& ub) {
  check_bounded<T, L, U>("lub_free", "Bounded variable", y, lb, ub);
  if (lb == NEGATIVE_INFTY) {
    return identity_constrain(ub_free(std::forward<T>(y), ub), lb);
  }
  if (ub == INFTY) {
    return identity_constrain(lb_free(std::forward<T>(y), lb), ub);
  }
  return logit((y - lb) / (ub - lb));
}

template <typename T, typename L, typename U, require_eigen_t<T>* = nullptr>
inline auto lub_free(T&& y, const L& lb, const U& ub) {
  check_bounded<T, L, U>("lub_free", "Bounded variable", y, lb, ub);
  if (lb == NEGATIVE_INFTY) {
    return identity_constrain(ub_free(std::forward<T>(y), ub), lb);
  }
  if (ub == INFTY) {
    return identity_constrain(lb_free(std::forward<T>(y), lb), ub);
  }
  return logit((y.array() - lb) / (ub - lb)).matrix().eval();
}

template <typename T, typename L, typename U,
          require_std_vector_t<T>* = nullptr>
inline auto lub_free(T&& y, const L& lb, const U& ub) {
  check_bounded<T, L, U>("lub_free", "Bounded variable", y, lb, ub);
  if (lb == NEGATIVE_INFTY) {
    return identity_constrain(ub_free(std::forward<T>(y), ub), lb);
  }
  if (ub == INFTY) {
    return identity_constrain(lb_free(std::forward<T>(y), lb), ub);
  }
  std::vector<return_type_t<T, L, U>> y_ret(y.size());
  const auto ub_lb_diff = ub - lb;
  for (int i = 0; i < y.size(); ++i) {
    y_ret[i] = logit((y[i] - lb) / ub_lb_diff);
  }
  return y_ret;
}

}  // namespace math
}  // namespace stan
#endif
