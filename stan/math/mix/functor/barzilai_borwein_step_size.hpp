#ifndef STAN_MATH_MIX_FUNCTOR_BARZILAI_BORWEIN_STEP_SIZE_HPP
#define STAN_MATH_MIX_FUNCTOR_BARZILAI_BORWEIN_STEP_SIZE_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace stan::math::internal {
/**
 * @brief  Curvature-aware Barzilai–Borwein (BB) step length with robust
 * safeguards.
 *
 * Given successive parameter displacements \f$s = x_k - x_{k-1}\f$ and
 * gradients \f$g_k\f$, \f$g_{k-1}\f$, this routine forms
 * \f$y = g_k - g_{k-1}\f$ and computes the two classical BB candidates
 *
 * \f{align*}{
 *   \alpha_{\text{BB1}} &= \frac{\langle s,s\rangle}{\langle s,y\rangle},\\
 *   \alpha_{\text{BB2}} &= \frac{\langle s,y\rangle}{\langle y,y\rangle},
 * \f}
 *
 * then chooses between them using the **spectral cosine**
 * \f$r = \cos^2\!\angle(s,y) = \dfrac{\langle s,y\rangle^2}
 *                                      {\langle s,s\rangle\,\langle
 * y,y\rangle}\in[0,1]\f$:
 *
 *  - if \f$r > 0.9\f$ (well-aligned curvature) and the previous line search
 *    did **≤ 1** backtrack, prefer the “long” step \f$\alpha_{\text{BB1}}\f$;
 *  - if \f$0.1 \le r \le 0.9\f$, take the neutral geometric mean
 *    \f$\sqrt{\alpha_{\text{BB1}}\alpha_{\text{BB2}}}\f$;
 *  - otherwise default to the “short” step \f$\alpha_{\text{BB2}}\f$.
 *
 * All candidates are clamped into \f$[\text{min\_alpha},\,\text{max\_alpha}]\f$
 * and must be finite and positive.
 * If the curvature scalars are ill-posed (non-finite or too small),
 * \f$\langle s,y\rangle \le \varepsilon\f$, or if `last_backtracks == 99`
 * (explicitly disabling BB for this iteration), the function falls back to a
 * **safe** step:
 * use `prev_step` when finite and positive, otherwise \f$1.0\f$, then clamp to
 * \f$[\text{min\_alpha},\,\text{max\_alpha}]\f$.
 *
 * @param s          Displacement between consecutive iterates
 *                   (\f$s = x_k - x_{k-1}\f$).
 * @param g_curr     Gradient at the current iterate \f$g_k\f$.
 * @param g_prev     Gradient at the previous iterate \f$g_{k-1}\f$.
 * @param prev_step  Previously accepted step length (used by the fallback).
 * @param last_backtracks
 *                   Number of backtracking contractions performed by the most
 *                   recent line search; set to 99 to force the safe fallback.
 * @param min_alpha  Lower bound for the returned step length.
 * @param max_alpha  Upper bound for the returned step length.
 *
 * @return A finite, positive BB-style step length \f$\alpha \in
 *         [\text{min\_alpha},\,\text{max\_alpha}]\f$ suitable for seeding a
 *         line search or as a spectral preconditioner scale.
 *
 * @note Uses \f$\varepsilon=10^{-16}\f$ to guard against division by very
 *       small curvature terms, and applies `std::abs` to BB ratios to avoid
 *       negative steps; descent is enforced by the line search.
 * @warning The vectors must have identical size. Non-finite inputs yield the
 *          safe fallback.
 */
inline double barzilai_borwein_step_size(const Eigen::VectorXd& s,
                                         const Eigen::VectorXd& g_curr,
                                         const Eigen::VectorXd& g_prev,
                                         double prev_step, int last_backtracks,
                                         double min_alpha, double max_alpha) {
  // Fallbacks
  auto safe_fallback = [&]() -> double {
    double a = std::clamp(
        prev_step > 0.0 && std::isfinite(prev_step) ? prev_step : 1.0,
        min_alpha, max_alpha);
    return a;
  };

  const Eigen::VectorXd y = g_curr - g_prev;
  const double sty = s.dot(y);
  const double sts = s.squaredNorm();
  const double yty = y.squaredNorm();

  // Basic validity checks
  constexpr double eps = 1e-16;
  if (!(std::isfinite(sty) && std::isfinite(sts) && std::isfinite(yty))
      || sts <= eps || yty <= eps || sty <= eps || last_backtracks == 99) {
    return safe_fallback();
  }

  // BB candidates
  double alpha_bb1 = std::clamp(std::abs(sts / sty), min_alpha, max_alpha);
  double alpha_bb2 = std::clamp(std::abs(sty / yty), min_alpha, max_alpha);

  // Safeguard candidates
  if (!std::isfinite(alpha_bb1) || !std::isfinite(alpha_bb2) || alpha_bb1 <= 0.0
      || alpha_bb2 <= 0.0) {
    return safe_fallback();
  }

  // Spectral cosine r = cos^2(angle(s, y)) in [0,1]
  const double r = (sty * sty) / (sts * yty);

  // Heuristic thresholds (robust defaults)
  constexpr double kLoose = 0.9;  // "nice" curvature
  constexpr double kTight = 0.1;  // "dodgy" curvature

  double alpha0 = alpha_bb2;  // default to short BB for robustness
  if (r > kLoose && last_backtracks <= 1) {
    // Spectrum looks friendly and line search was not harsh -> try long BB
    alpha0 = alpha_bb1;
  } else if (r >= kTight && r <= kLoose) {
    // Neither clearly friendly nor clearly dodgy -> neutral middle
    alpha0 = std::sqrt(alpha_bb1 * alpha_bb2);
  }  // else keep alpha_bb2

  // Clip to user bounds
  alpha0 = std::clamp(alpha0, min_alpha, max_alpha);

  if (!std::isfinite(alpha0) || alpha0 <= 0.0) {
    return safe_fallback();
  }
  return alpha0;
}

}  // namespace stan::math::internal
#endif
