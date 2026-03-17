#ifndef STAN_MATH_PRIM_CONSTRAINT_CORR_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_CORR_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/atanh.hpp>

namespace stan {
namespace math {

/**
 * Return the unconstrained scalar that when transformed to
 * a valid correlation produces the specified value.
 *
 * <p>This function inverts the transform defined for
 * <code>corr_constrain(T)</code>, which is the inverse hyperbolic
 * tangent,
 *
 * <p>\f$ f^{-1}(y)
 *          = \mbox{atanh}\, y
 *          = \frac{1}{2} \log \frac{y + 1}{y - 1}\f$.
 *
 * @tparam T Type of correlation
 * @param[in] y correlation
 * @return free scalar that transforms to the specified input
 */
template <typename T>
inline plain_type_t<T> corr_free(T&& y) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  check_bounded("lub_free", "Correlation variable", y_ref, -1.0, 1.0);
  return atanh(std::forward<decltype(y_ref)>(y_ref));
}

}  // namespace math
}  // namespace stan
#endif
