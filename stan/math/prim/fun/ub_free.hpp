#ifndef STAN_MATH_PRIM_FUN_UB_FREE_HPP
#define STAN_MATH_PRIM_FUN_UB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/log.hpp>
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
 * @tparam T type of scalar or matrix
 * @tparam U type of upper bound
 * @param y constrained scalar with specified upper bound
 * @param ub upper bound
 * @return unconstrained scalar with respect to upper bound
 * @throw std::invalid_argument if constrained scalar is greater
 *   than the upper bound.
 */
template <typename T, typename U>
inline auto ub_free(T&& y, U&& ub) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  auto&& ub_ref = to_ref(std::forward<U>(ub));
  using y_t = decltype(y_ref);
  using ub_t = decltype(ub_ref);
  check_less_or_equal("ub_free", "Upper bounded variable", value_of(y_ref),
                      value_of(ub_ref));
  if (is_positive_infinity(ub_ref)) {
    return identity_constrain(std::forward<y_t>(y_ref), ub_ref);
  } else {
    return log(subtract(std::forward<ub_t>(ub_ref), std::forward<y_t>(y_ref)));
  }
}

}  // namespace math
}  // namespace stan
#endif
