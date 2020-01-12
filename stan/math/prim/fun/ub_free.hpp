#ifndef STAN_MATH_PRIM_FUN_UB_FREE_HPP
#define STAN_MATH_PRIM_FUN_UB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the free scalar that corresponds to the specified
 * upper-bounded value with respect to the specified upper bound.
 *
 * <p>The transform is the reverse of the
 * <code>ub_constrain(T, double)</code> transform,
 *
 * <p>\f$f^{-1}(y) = \log -(y - U)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * If the upper bound is positive infinity, this function
 * reduces to <code>identity_free(y)</code>.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @param y constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename T, typename U>
inline return_type_t<T, U> ub_free(const T& y, const U& ub) {
  using std::log;
  if (ub == INFTY) {
    return identity_free(y);
  }
  check_less_or_equal("ub_free", "Upper bounded variable", y, ub);
  return log(ub - y);
}

}  // namespace math
}  // namespace stan
#endif
