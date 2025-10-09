#ifndef STAN_MATH_MIX_FUNCTOR_WOLFE_LINE_SEARCH_HPP
#define STAN_MATH_MIX_FUNCTOR_WOLFE_LINE_SEARCH_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/wolfe_line_search.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
// #define LAPLACE_DEBUG
#ifdef LAPLACE_DEBUG
#include <iomanip>
#include <iostream>
#include <ostream>
#include <optional>
#include <fstream>
#endif



namespace stan::math {

  /**
 * @brief Options for Wolfe line search during optimization.
 *
 * These parameters control the behavior of the line search,
 * including the Wolfe condition thresholds, step size bounds,
 * and numerical tolerances.
 */
struct laplace_line_search_options {
  /**
   * @brief Armijo condition parameter (sufficient decrease).
   *
   * Controls how much the objective function must decrease
   * relative to the step length and gradient. Must be in (0, 1).
   * Smaller values make the condition easier to satisfy.
   */
  double c1{1e-4};

  /**
   * @brief Curvature condition parameter.
   *
   * Ensures that the step length provides a sufficient reduction
   * in the gradient magnitude. Must satisfy c1 < c2 < 1.
   * Larger values enforce stronger curvature requirements.
   */
  double c2{0.9};

  /**
   * @brief Backtracking shrinkage factor.
   *
   * When the Wolfe conditions are not met, the step size is
   * multiplied by this factor (0 < tau < 1) to reduce it.
   */
  double tau{0.5};

  /**
   * @brief Minimum allowable step size.
   *
   * Prevents the line search from shrinking the step length
   * below this threshold, which helps avoid numerical issues.
   */
  double min_alpha{1e-8};

  /**
   * @brief Maximum allowable step size.
   *
   * Caps the growth of the step length to ensure stability.
   */
  double max_alpha{32.0};

  /**
   * @brief Step size expansion factor.
   *
   * If the Wolfe conditions are satisfied and further progress
   * seems possible, the step size may be multiplied by this factor
   * to accelerate convergence.
   */
  double scale_up{2.0};

  /**
   * @brief Absolute gradient tolerance.
   *
   * If the gradient norm falls below this threshold, the line
   * search may terminate early since further progress is negligible.
   */
  double abs_grad_threshold{1e-12};

  /**
   * @brief Absolute objective tolerance.
   *
   * If the change in objective function value falls below this
   * threshold, the line search may terminate as improvements are
   * considered insignificant.
   */
  double abs_obj_threshold{1e-12};
};
namespace internal {
namespace debug {

constexpr void print(int tabs) {}

#ifdef LAPLACE_DEBUG
template <typename Val, typename... Types>
void print(int tabs, const char* name, Val&& val, Types&&... args) {
  csv_file << name << ", " << std::setprecision(13) << std::fixed
           << stan::math::eval(val) << "\n";
  print(tabs, std::forward<Types>(args)...);
}

template <typename Val, typename... Types>
void print(const char* msg, int tabs, const char* name, Val&& val,
           Types&&... args) {
  csv_file << name << ", " << std::setprecision(13) << std::fixed
           << stan::math::eval(val) << "\n";
  print(tabs, std::forward<Types>(args)...);
}

void print(const char* msg) {
  // std::cout << msg << "\n";
}

template <typename Val>
void print(const char* msg, Val&& val) {
  csv_file << msg << ", " << std::setprecision(13) << std::fixed
           << stan::math::eval(val) << "\n";
}

#else
template <typename Val, typename... Types>
constexpr void print(int tabs, const char* name, Val&& val, Types&&... args) {}

template <typename Val, typename... Types>
constexpr void print(const char* msg, int tabs, const char* name, Val&& val,
                     Types&&... args) {}

constexpr void print(const char* msg) {}

template <typename Val>
constexpr void print(const char* msg, Val&& val) {}
#endif
}  // namespace debug





/*
* Safeguarded cubic-or-bisection step chooser for MAXIMIZATION.
* Given a bracket [a, b] with values and directional derivatives
* (fa, fpa) at a and (fb, fpb) at b, returns a trial point in (a, b).
*
* Preconditions (recommended for a well-formed maximum bracket):
*   a < b, fpa > 0, fpb < 0.
* If inputs are degenerate or non-finite, falls back to the midpoint.
*
* The routine internally normalizes the interval to s in [0, 1], fits a
* cubic Hermite model F(s), and selects among:
*   - stationary points of F(s) (cubic argmax),
*   - a secant estimate for the derivative root,
*   - pure bisection (midpoint),
* then clamps away from edges by `edge_guard`.
*/
template <typename Scalar>
inline Scalar cubic_or_bisect_max(Scalar a, Scalar fa, Scalar fpa,
                           Scalar b, Scalar fb, Scalar fpb) {
  // Basic validation.
  if (!(b > a) || !std::isfinite(fa) || !std::isfinite(fb) ||
      !std::isfinite(fpa) || !std::isfinite(fpb)) {
    return (a + b) / Scalar(2.0);
  }

  const Scalar width = b - a;
  const Scalar fpa_s = width * fpa;  // slope w.r.t. s at s=0
  const Scalar fpb_s = width * fpb;  // slope w.r.t. s at s=1

  // Cubic Hermite coefficients in s \in [0,1]:
  // F(s) = a3*s^3 + a2*s^2 + a1*s + a0
  // with F(0)=fa, F'(0)=fpa_s, F(1)=fb, F'(1)=fpb_s.
  const Scalar a0 = fa;
  const Scalar a1 = fpa_s;
  const Scalar a2 = (Scalar(3) * (fb - fa)) - (Scalar(2) * fpa_s) - fpb_s;
  const Scalar a3 = (Scalar(2) * (fa - fb)) + fpa_s + fpb_s;

  auto eval = [&](Scalar s) -> Scalar {
    // Horner evaluation of F(s).
    return ((a3 * s + a2) * s + a1) * s + a0;
  };

  // Helper: push candidate s if it's inside (0,1).
  struct Candidate { Scalar s; Scalar val; };
  Candidate best{Scalar(0.5), eval(Scalar(0.5))};  // start with bisection

  auto consider = [&](Scalar s) {
    if (!(s > Scalar(0) && s < Scalar(1)) || !std::isfinite(s)) return;
    const Scalar v = eval(s);
    if (v > best.val) best = {s, v};
  };

  // 1) Stationary points of the cubic model (argmax/min candidates).
  // F'(s) = 3*a3*s^2 + 2*a2*s + a1 = 0.
  {
    const Scalar A = Scalar(3) * a3;
    const Scalar B = Scalar(2) * a2;
    const Scalar C = a1;

    const Scalar eps = std::numeric_limits<Scalar>::epsilon();
    if (std::abs(A) <= eps * (std::abs(B) + std::abs(C))) {
      // Degenerates to linear: B*s + C = 0.
      if (std::abs(B) > eps * std::abs(C)) {
        consider(-C / B);
      }
    } else {
      const Scalar disc = std::fma(-Scalar(4) * A, C, B * B);  // B^2 - 4AC
      if (disc >= 0) {
        const Scalar r = std::sqrt(disc);
        // q-formula for numerical stability.
        const Scalar q = -Scalar(0.5) * (B + std::copysign(r, B));
        const Scalar s1 = q / A;
        const Scalar s2 = (q == Scalar(0)) ? std::numeric_limits<Scalar>::quiet_NaN()
                                           : C / q;
        consider(s1);
        consider(s2);
      }
    }
  }

  // 2) Secant estimate for the derivative root (More–Thuente often uses
  //    a secant/minimizer mix; we add the secant root for robustness).
  // Derivative across s is linearly interpolated from fpa to fpb:
  // D(s) ~ fpa + (fpb - fpa)*s; root at s = fpa / (fpa - fpb).
  if (std::isfinite(fpa) && std::isfinite(fpb) && (fpa != fpb)) {
    const Scalar s_secant = fpa / (fpa - fpb);
    consider(s_secant);
  }

  // Edge guard: keep away from exact ends.
  constexpr Scalar edge_guard = Scalar(1e-3);
  constexpr Scalar lo = 0.0 * (1.0 - edge_guard) + (1.0 * edge_guard);
  constexpr Scalar hi = Scalar(0.0) * (1.0 - (Scalar(1) - edge_guard)) + (1.0 * (Scalar(1) - edge_guard));
  const Scalar s_star = std::clamp(best.s, lo, hi);
  // Map back to alpha-space.
  return a + s_star * width;
}

template <typename Eval, typename Options>
inline auto cubic_or_bisect_max(Eval&& low, Eval&& high, Options&& opt) {
  auto alpha = cubic_or_bisect_max(low.alpha, low.obj, low.dir, high.alpha, high.obj, high.dir);
  alpha
        = std::clamp(alpha, opt.line_search.min_alpha, high.alpha * 0.9);
  return alpha;

}

template <typename Option>
inline auto check_armijo(double obj_next, double obj_init, double alpha_next,
                         double dir0, Option&& opt) {
  debug::print(
      "check_armijo: ", 2, "armijo:    ",
      (obj_next >= obj_init + alpha_next * dir0 * opt.line_search.c1 ? 1 : 0),
      "obj_next:   ", obj_next, "obj_init:   ", obj_init,
      "alpha_next: ", alpha_next, "dir0:       ", dir0,
      "c1:         ", opt.line_search.c1, "obj + alpha * dir0 * c1: ",
      (obj_init + alpha_next * dir0 * opt.line_search.c1));
  return obj_next >= obj_init + alpha_next * dir0 * opt.line_search.c1;
}

template <typename Option>
inline auto check_wolfe_curve(double dir_deriv_next, double dir_deriv_init,
                              Option&& opt) {
  debug::print("check_wolfe_curve: ", 2, "wolfe:    ",
               (std::abs(dir_deriv_next)
                        <= (opt.line_search.c2 * std::abs(dir_deriv_init))
                    ? 1
                    : 0),
               "deriv_next: ", dir_deriv_next, "deriv_init: ", dir_deriv_init,
               "c2:         ", opt.line_search.c2,
               "abs(d_next):   ", std::abs(dir_deriv_next), "abs(d_init)*c2 ",
               (std::abs(dir_deriv_init) * opt.line_search.c2));
  return std::abs(dir_deriv_next)
         <= (opt.line_search.c2 * std::abs(dir_deriv_init));
}



enum class wolfe_return {
  PASS,         // Success
  ARMIJO_PASS,  // Armijo condition passed but Wolfe failed
  FAIL,         // Armijo Wolfe conditions failed
  GRAD_CONV     // Approximate gradient converged
};

/**
 * @brief  Strong‑Wolfe line search with cubic‑interpolation "zoom" for
 *Laplace‐style log‑likelihood problems.
 *
 * This routine searches along the space of the latent gaussian *a*
 * \f$a(\alpha) = a_prev + \alpha p`, `p = a − a_prev\f$,
 * looking for the largest step \f$\alpha\f$ that satisfies the **strong‑Wolfe**
 * conditions
 *
 * \f{align*}{
 *   \phi(\alpha) &\;\ge\; \phi(0) + c_1 \alpha \phi'(0) \quad\text{(Armijo)},\\
 *   |\phi'(\alpha)| &\;\le\; c_2 |\phi'(0)| \quad\text{(curvature)}, \f}
 *
 * where
 *  \f$\phi(\alpha)=\text{obj\_fun}\bigl(a(\alpha),\;\theta(\alpha)\bigr)\f$,
 *  \f$\theta(\alpha)=\text{covariance}\;·\;a(\alpha)\f$,
 *  \f$\phi'(\alpha)=\nabla\phi(\alpha)^{\!T} p\f$.
 *
 * The search proceeds in three phases
 *
 *  1. **Back‑tracking** – halve the initial \f$\alpha\f$ until both Wolfe
 *conditions pass.
 *  2. If the \f$alpha\f$ is 1 or goes below `min_alpha`, then we end the search
 *early, as Laplace problems commonly accept a full Newton step.
 *  3. **Bracketing by doubling** – starting from that good \f$\alpha\f$, double
 *the step until Armijo fails; the last good point is the left end of the
 *bracket, the first failing point the right end.
 *  4. **Cubic zoom** – repeatedly fit a cubic through the bracket end‑points
 *     (`cubic_interp`), evaluate the objective/gradient at the
 *     predicted maximiser, and shrink the bracket until a Wolfe‑compliant
 *     step is found or the interval width falls below
 *`opt.line_search.min_alpha`.
 *
 * * **Gradient reuse** – the caller provides `theta_grad` computed at the
 *   starting point; inside the loop it is overwritten with fresh gradients
 *   via `laplace_likelihood::theta_grad`.
 * * **Early‑exit for \f$\alpha\f$ = 1 or < `min_alpha`** Laplace problems
 *commonly accept a full Newton step. The function short‑circuits to avoid any
 *extra work.
 *
 * @tparam F Callable type of the raw log‑likelihood (passed to
 *`laplace_likelihood::theta_grad`)
 * @tparam Obj Callable returning the scalar objective value `obj_fun(a, theta)`
 * @tparam Grad Callable returning the gradient in *a‑space* given `(a, theta,
 *theta_grad)`
 * @tparam LLArgs Struct or tuple holding additional arguments forwarded to
 *`ll_fun` and `theta_grad`
 * @tparam Stream Any type that implements the stream interface used by
 *`laplace_likelihood` for diagnostic messages; may be `std::ostream`‐like or
 *`nullptr`
 * @tparam Options Struct holding search parameters (`c1`, `c2`, `min_alpha`,
 *...)
 *
 * @param[in,out] theta Current \f$\theta=\Sigma a\f$; overwritten with the
 *accepted value on success
 * @param[in,out] obj_init Objective at \f$\alpha\f$ = 0 on entry; updated to
 *the objective at the accepted step on exit
 * @param[in,out] alpha_init Step size suggestion on input; on success holds the
 *step that satisfied strong‑Wolfe
 * @param[in,out] theta_grad Gradient of the log‑likelihood wrt theta at the
 ***current** theta.  Recomputed internally and returned at the final theta
 * @param[in,out] a Working copy of *a*; on success contains the accepted point
 * @param[in] a_prev Starting point (\f$\alpha\f$ = 0)
 * @param[in] ll_fun Raw log‑likelihood functor
 * @param[in] obj_fun Objective functor to be maximised
 * @param[in] grad_fun Functor returning gradint of objective with respect to
 *`a`.
 * @param[in] covariance Symmetric positive‑definite $\Sigma$ converting *a* to
 *theta
 * @param[in] llzfif _args Extra arguments forwarded to `ll_fun` and
 *`theta_grad`
 * @param[in] opt Line‑search constants (`c1`, `c2`, etc.)
 * @param[in,out] msgs Optional diagnostics stream.  May be `nullptr` to
 *suppress messages
 *
 * @return `true`  if a step satisfying both Wolfe conditions was found
 *         `false` if only Armijo is satisfied (strong‑Wolfe failed but a
 *                “least‑bad’’ step was returned).
 *
 * @warning The helper `cubic_interp` assumes its first point
 *          corresponds to *g*(0)=0; do **not** alter the baseline
 *          initialisation logic or the interpolation will become
 *          inconsistent.
 */
template <typename F, class Obj, class Grad, typename LLArgs, typename Stream,
          typename Options>
inline wolfe_return wolfe_line_search(
    Eigen::VectorXd& theta, double& obj_init, double& alpha_init,
    Eigen::VectorXd& theta_grad, Eigen::VectorXd& a,
    const Eigen::VectorXd& a_prev, F&& ll_fun, Obj&& obj_fun, Grad&& grad_fun,
    const Eigen::MatrixXd& covariance, LLArgs&& ll_args, Options&& opt,
    Stream* msgs) {
  struct Eval {
    double alpha;   // alpha
    double obj;   // obj
    double dir;   // directional derivative
  };
  Eigen::VectorXd p = a - a_prev;
  double dir_deriv_init = grad_fun(a_prev, theta, theta_grad).dot(p);
  Eval low{0.0, obj_init, dir_deriv_init};
  Eval best = low;  // keep the best Armijo-OK in case strong-Wolfe fails

  auto armijo_ok = [&](const Eval& eval) -> bool {
    return check_armijo(eval.obj, obj_init, eval.alpha, dir_deriv_init, opt);
  };
  auto wolfe_ok = [&](const Eval& eval) -> bool {
    return check_wolfe_curve(eval.dir, dir_deriv_init, opt);
  };
  int total_updates = 0;
  auto update_step
      = [&p, &a_prev, &covariance, &ll_fun, &ll_args, &obj_fun, &grad_fun, msgs,
         &total_updates](auto& a_in, auto& theta_in, auto& theta_grad_in,
                         auto& eval_in) {
          total_updates++;
          a_in = a_prev + eval_in.alpha * p;
          theta_in = covariance * a_in;
          theta_grad_in
              = laplace_likelihood::theta_grad(ll_fun, theta_in, ll_args, msgs);
          eval_in.obj = obj_fun(a_in, theta_in);
          eval_in.dir = grad_fun(a_in, theta_in, theta_grad_in).dot(p);
        };
  double alpha_start = std::clamp(alpha_init * 2.0, opt.line_search.min_alpha,
                                 opt.line_search.max_alpha);
  Eigen::VectorXd a_try(a_prev.size());
  Eigen::VectorXd theta_try(a_prev.size());
  Eval high{alpha_start, obj_init, dir_deriv_init};
  update_step(a_try, theta_try, theta_grad, high);
  {
    while (!(std::isfinite(high.obj) && theta_try.allFinite())) {
      high.alpha *= 0.5;
      if (high.alpha < opt.line_search.min_alpha) {
        debug::print("Exit on precheck numerical trouble", 1);
        debug::print("total_updates", total_updates);
        alpha_init = high.alpha;
        return wolfe_return::FAIL;
      }
      update_step(a_try, theta_try, theta_grad, high);
    }
    debug::print("First precheck: ", 1, "high.alpha: ", high.alpha,
                 "high.obj:   ", high.obj, "deriv_high: ", high.dir,
                 "deriv_init: ", dir_deriv_init);
    // Quick accept if Armijo and Wolfe conditions are satisfied
    if (armijo_ok(high)) {
      if (wolfe_ok(high)) {
        debug::print("Exit on first precheck", 1);
        a.swap(a_try);
        theta.swap(theta_try);
        obj_init = high.obj;
        alpha_init = high.alpha;
        debug::print("total_updates", total_updates);
        return wolfe_return::PASS;
      } else {
        if (best.obj < high.obj) {
          best = high;
        }
      }
    }
  }
  // If current alpha fails, backtrack down till we find a good point
  debug::print("Begin Loop: ", 1, "Initial alpha: ", high.alpha);
  //    "g0:            ", g0.transpose().eval(),
  //    "theta_try:     ", theta_try.transpose().eval());
  int loop_iter = 0;
  const auto grad_tol = opt.line_search.abs_grad_threshold;
  const auto obj_tol = opt.line_search.abs_obj_threshold;
  // If true we have already found a good first point
  bool found_right = false;

  /**
   * For each case
   * armijo | wolfe | sign(g) | Action
   * -------+-------+---------+---------------------------------------------
   * [1]  T     |   T   |         | Accept α (PASS: strong-Wolfe satisfied)
   * [2]  T     |   F   |   > 0   | Promote left: set α_low ← α, keep best, expand α (e.g. α ← scale_up·α)
   * [3]  T     |   F   |   < 0   | Set α_high ← α, stop expanding → zoom
   * [4]  F     |   T   |         | Set α_high ← α, stop expanding → zoom
   * [5]  F     |   F   |         | Set α_high ← α, stop expanding → zoom
   **/
  while (!found_right && high.alpha < opt.line_search.max_alpha) {
    // 1. Evaluate f(α) and g(α)
    update_step(a_try, theta_try, theta_grad, high);
    debug::print("First While", 1, "Second Iter:       ", loop_iter++,
                 "high.alpha: ", high.alpha, "high.obj:   ", high.obj,
                 "deriv_high: ", high.dir, "deriv_init: ", dir_deriv_init,
                 "theta_try:  ", theta_try.transpose());
    const bool finite_ok = std::isfinite(high.obj) && theta_try.allFinite();
    // 2. Handle numerical trouble first
    if (!finite_ok) {  //   f or g is NaN/Inf → shrink
      high.alpha *= 0.5;
      if (high.alpha < opt.line_search.min_alpha) {
        break;
      }
      continue;
    }
    if (armijo_ok(high)) {
      // [1]
      if (wolfe_ok(high)) {
        a.swap(a_try);
        theta.swap(theta_try);
        obj_init = high.obj;
        alpha_init = high.alpha;
        debug::print("Exit on first while", 1);
        debug::print("total_updates", total_updates);
        return wolfe_return::PASS;
      } else {
        if (best.obj < high.obj) {
          best = high;
        }
        // [2]
        if (std::signbit(high.dir)) {
          low = high;
          high.alpha *= opt.line_search.scale_up;
          continue;
        } else {
          // [3]
          found_right = true;
        }
      }
    }
    // [4,5]
    found_right = true;
    break;
  }
  // Check for grad convergence
  if (std::abs(high.dir) <= grad_tol ||  // tiny slope
    std::abs(high.obj - obj_init) <= obj_tol
        && high.alpha < 1e-4) {                // tiny gain
    if (std::abs(high.dir) <= grad_tol &&  // tiny slope
        std::abs(high.obj - obj_init) <= obj_tol) {
      debug::print("Exit on grad_tol and obj_tol", 1);
    } else if (std::abs(high.dir) <= grad_tol) {
      debug::print("Exit on grad_tol", 1);
    } else if (std::abs(high.obj - obj_init) <= obj_tol) {
      debug::print("Exit on obj_tol", 1);
    } else {
      debug::print("Exit on alpha failure", 1);
    }
    a.swap(a_try);
    theta.swap(theta_try);
    obj_init = high.obj;
    alpha_init = high.alpha;
    debug::print("total_updates", total_updates);
    return wolfe_return::GRAD_CONV;  // add a code in your enum
  }

  update_step(a_try, theta_try, theta_grad, high);

  debug::print("_______End First While: ", 1, "high.alpha: ", high.alpha,
               "high.obj:   ", high.obj, "high.dir: ", high.dir,
               "dir_deriv_init: ", dir_deriv_init);
  loop_iter = 0;
  Eval mid{low};
  while (high.alpha - low.alpha > opt.line_search.min_alpha) {
    const double diff_alpha = high.alpha - low.alpha;
    mid.alpha = cubic_or_bisect_max(low, high, opt);
    update_step(a_try, theta_try, theta_grad, mid);
    debug::print("Cube: ", 1, "Cube Iter:           ", loop_iter++,
                 "mid.alpha:      ", mid.alpha, "mid.obj:        ", mid.obj,
                 "high.dir: ", mid.dir,
                 "low.alpha:      ", low.alpha, "low.obj:        ", low.obj,
                 "low.dir:  ", low.dir,
                 "high.alpha:     ", high.alpha, "high.obj:       ", high.obj,
                 "high.dir: ", high.dir);
    const bool finite_ok = std::isfinite(mid.obj) && theta_try.allFinite();
    if (!finite_ok) {
      debug::print("Exit on failed finite test", 1);
      high = mid;
      continue;
    }
    /**
     * | Armijo | Wolfe | mid.obj > low.obj | mid.dir * low.dir < 0 | Action |
     * | T     |   T   |                   |                       | Accept |
     * | T     |   F   |         T         |           T           | Dir Sign change pins high = mid |
     * | T     |   F   |         T         |           F           | low = mid |
     * | T     |   F   |         F         |           F           | low = mid |
     * | T     |   F   |         F         |           F           | high = mid |
     * | F     |   F   |         F         |           F           | high = mid |
     */
    if (armijo_ok(mid)) {
      if (wolfe_ok(mid)) {
        a.swap(a_try);
        theta.swap(theta_try);
        obj_init = mid.obj;
        alpha_init = mid.alpha;
        debug::print("Exit on safe on zoom", 1);
        debug::print("total_updates", total_updates);
        return wolfe_return::PASS;
      } else {
        if (best.obj < mid.obj) {
          best = mid;
        }
        if (mid.obj > low.obj) {
          // sign change
          if (mid.dir * low.dir < 0) {
            high = mid;
          } else {
            low = mid;
          }
        } else {
          high = mid;
        }
      }
    } else {
      high = mid;
    }
    if (std::abs(mid.dir) <= grad_tol
        || std::abs(mid.obj - obj_init) <= obj_tol && mid.alpha < 1e-6) {
      if (std::abs(mid.dir) <= grad_tol &&  // tiny slope
          std::abs(mid.obj - obj_init) <= obj_tol) {
        debug::print("Exit on grad_tol and obj_tol", 1);
      } else if (std::abs(mid.dir) <= grad_tol) {
        debug::print("Exit on grad_tol", 1);
      } else if (std::abs(mid.obj - obj_init) <= obj_tol) {
        debug::print("Exit on obj_tol", 1);
      } else {
        debug::print("Exit on cube alpha failure", 1);
      }
      a.swap(a_try);
      theta.swap(theta_try);
      obj_init = mid.obj;
      alpha_init = mid.alpha;
      debug::print("total_updates", total_updates);
      return wolfe_return::GRAD_CONV;
    }
  }
  debug::print("Failed zoom: ", 1, "Failed zoom:", 1, "low.alpha: ", low.alpha,
               "low.obj:   ", low.obj, "deriv_low: ", low.dir,
               "high.alpha:", high.alpha, "high.obj:  ", high.obj,
               "deriv_high:", high.dir);
  // On failure, use the best point we have found so far that at least satisfies armijo
  update_step(a_try, theta_try, theta_grad, best);
  a.swap(a_try);
  theta.swap(theta_try);
  theta_grad = laplace_likelihood::theta_grad(ll_fun, theta_try, ll_args, msgs);
  // We already calculated best.obj and theta_grad so no need to recompute here
  best.dir = grad_fun(a, theta, theta_grad).dot(p);
  const bool armijo_ok_mid
      = check_armijo(best.obj, obj_init, best.alpha, dir_deriv_init, opt);
  const bool curve_ok_mid = check_wolfe_curve(best.dir, dir_deriv_init, opt);
  obj_init = best.obj;
  alpha_init = best.alpha;
  if (armijo_ok_mid && curve_ok_mid) {
    debug::print("Exit on safe after zoom", 1);
    debug::print("total_updates", total_updates);
    return wolfe_return::PASS;
  } else if (armijo_ok_mid) {
    debug::print("Exit on only satisfying armijo", 1);
    debug::print("total_updates", total_updates);
    return wolfe_return::ARMIJO_PASS;
  } else {
    debug::print("Exit on failure", 1);
    debug::print("total_updates", total_updates);
    return wolfe_return::FAIL;
  }
}
}

}
#endif
