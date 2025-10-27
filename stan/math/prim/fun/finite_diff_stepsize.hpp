#ifndef STAN_MATH_PRIM_FUN_FINITE_DIFF_STEPSIZE_HPP
#define STAN_MATH_PRIM_FUN_FINITE_DIFF_STEPSIZE_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <cmath>

namespace stan {
namespace math {

#include <cmath>
#include <limits>
#include <algorithm>

namespace internal {
template <std::size_t StencilOrder = 2, typename T>
inline constexpr auto eps_root_calc() {
  constexpr T eps = std::numeric_limits<T>::epsilon();
  if constexpr (StencilOrder == 2) {
    // ε^(1/3)
    return std::cbrt(eps);
  } else if constexpr (StencilOrder == 3) {
    // ε^(1/4) = sqrt(sqrt(ε))
    return std::sqrt(std::sqrt(eps));
  } else if constexpr (StencilOrder == 5) {
    // ε^(1/6) = sqrt(cbrt(ε))
    return std::sqrt(std::cbrt(eps));
  } else {
    // General fallback: ε^(1/(p+1))
    return std::pow(eps, T(1) / T(StencilOrder + 1));
  }
}
}
/**
 * @brief Compute a finite-difference step size suitable for a stencil with
 *        leading truncation order \p StencilOrder.
 *
 * This implements the standard balance between truncation error
 * \f$A\,h^{p}\f$ and floating-point rounding error \f$B\,\varepsilon/h\f$,
 * whose minimizer is \f$h_\star \propto \varepsilon^{1/(p+1)}\f$, where
 * \f$p = \texttt{StencilOrder}\f$ and \f$\varepsilon\f$ is machine epsilon
 * for the scalar type \p T. We additionally scale by \f$\max(1, |u|)\f$ to
 * obtain an absolute step size near the point being perturbed.
 *
 * For the common second-order (3-point central) stencil, the function uses an
 * `if constexpr` branch to compute \f$\varepsilon^{1/3}\f$ via `std::cbrt`
 * (avoids a `pow` and is numerically tidy). For higher orders it uses
 * \f$\varepsilon^{1/(p+1)}\f` via `std::pow`.
 *
 * @tparam T            Floating-point scalar type (e.g., float, double).
 * @tparam StencilOrder Leading truncation order \f$p\f$ of your stencil:
 *                      - first-derivative 3-point central  → \f$p=2\f$
 *                      - first-derivative 5-point central  → \f$p=4\f$
 *                      - first-derivative 6-point central  → \f$p=6\f$
 *                      - “derivative of Hessian” via 5-point diff → \f$p=4\f$
 *
 * @param u  The coordinate value (or local scale) at which the step will be
 *           applied; the step is scaled by \f$\max(1, |u|)\f$.
 * @param c  Dimensionless tuning constant multiplying the theoretically
 *           optimal step. Default is 1.0. In practice, \f$c \in [0.5, 2]\f$
 *           works well:
 *           - Increase \p c (larger h) if round-off/cancellation dominates
 *             (noisy function values, very large |f|, many subtractions).
 *           - Decrease \p c (smaller h) if truncation dominates (very smooth
 *             function with large higher derivatives or low curvature scale).
 *           If your problem has a known physical scale S (not |u|), prefer
 *           passing \p u scaled by S or modify the scaling accordingly.
 *
 * @return A step size \f$h\f$ such that \f$u+h \neq u\f$; if the theoretical
 *         step underflows at \p u, the function falls back to the next
 *         representable increment.
 *
 * @note The step computed here assumes smooth (at least \(p{+}1\) times
 *       differentiable) behavior along the perturbed coordinate.
 */
template <std::size_t StencilOrder = 2, typename T>
inline T finite_diff_stepsize(T u, T c = T(1)) {
  const T scale = std::max(T(1), std::abs(u));
  const T eps_root = internal::eps_root_calc<StencilOrder, T>();
  const T h = c * eps_root * scale;
  // Ensure perturbation isn’t rounded away at u.
  if (u + h == u) {
    const T next = std::nextafter(u, std::numeric_limits<T>::infinity());
    return std::max(h, next - u);
  }
  return h;
}

}  // namespace math
}  // namespace stan
#endif
