#ifndef STAN_MATH_PRIM_CONSTRAINT_PROB_FREE_HPP
#define STAN_MATH_PRIM_CONSTRAINT_PROB_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/logit.hpp>

namespace stan {
namespace math {

/**
 * Return the free scalar that when transformed to a probability
 * produces the specified scalar.
 *
 * <p>The function that reverses the constraining transform
 * specified in <code>prob_constrain(T)</code> is the logit
 * function,
 *
 * <p>\f$f^{-1}(y) = \mbox{logit}(y) = \frac{1 - y}{y}\f$.
 *
 * @tparam T type of constrained value
 * @param y constrained value
 * @return corresponding unconstrained value
 * @throw std::domain_error if y is not in (0, 1)
 */
template <typename T>
inline auto prob_free(T&& y) {
  auto&& y_ref = to_ref(std::forward<T>(y));
  check_bounded<plain_type_t<T>, double, double>(
      "prob_free", "Probability variable", y_ref, 0, 1);
  return logit(std::forward<decltype(y_ref)>(y_ref));
}

}  // namespace math
}  // namespace stan
#endif
