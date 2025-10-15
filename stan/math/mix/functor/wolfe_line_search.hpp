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
  double max_alpha{16.0};

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
        = std::clamp(alpha, opt.min_alpha, high.alpha * 0.9);
  return alpha;

}

template <typename Option>
inline auto check_armijo(double obj_next, double obj_init, double alpha_next,
                         double dir0, Option&& opt) {
  debug::print(
      "check_armijo: ", 2, "armijo:    ",
      (obj_next >= obj_init + alpha_next * dir0 * opt.c1 ? 1 : 0),
      "obj_next:   ", obj_next, "obj_init:   ", obj_init,
      "alpha_next: ", alpha_next, "dir0:       ", dir0,
      "c1:         ", opt.c1, "obj + alpha * dir0 * c1: ",
      (obj_init + alpha_next * dir0 * opt.c1));
  return obj_next >= obj_init + alpha_next * dir0 * opt.c1;
}

template <typename Option>
inline auto check_wolfe_curve(double dir_deriv_next, double dir_deriv_init,
                              Option&& opt) {
  debug::print("check_wolfe_curve: ", 2, "wolfe:    ",
               (std::abs(dir_deriv_next)
                        <= (opt.c2 * std::abs(dir_deriv_init))
                    ? 1
                    : 0),
               "deriv_next: ", dir_deriv_next, "deriv_init: ", dir_deriv_init,
               "c2:         ", opt.c2,
               "abs(d_next):   ", std::abs(dir_deriv_next), "abs(d_init)*c2 ",
               (std::abs(dir_deriv_init) * opt.c2));
  return std::abs(dir_deriv_next)
         <= (opt.c2 * std::abs(dir_deriv_init));
}

enum class WolfeReturn : uint8_t {
  // both conditions true
  Wolfe,
  // Armijo true, curvature false
  Armijo,
  // |phi'| small
  ConvergedGradient,
  // |phi - phi0| small
  ConvergedObjective,
  // |phi - phi0| and |phi'| small
  ConvergedObjectiveAndGradient,
  // high and low became too small
  IntervalTooSmall,
  // alpha < min_alpha
  StepTooSmall,
  // max iters reached
  ReachedMaxStep,
  // failed to find a step
  NumericalIssue,
  // All other failures
  Fail
};

struct WolfeStatus {
  // total updates/evaluations
  int num_evals_{0};
  int num_backtracks_{0};
  WolfeReturn stop_{WolfeReturn::Fail};
  WolfeStatus() = default;
  WolfeStatus(WolfeReturn stop, int evals, int back) :
    stop_(stop), num_evals_(evals), num_backtracks_(back) {}
};

inline auto wolfe_status_str(WolfeStatus s) {
  switch (s.stop_) {
  case WolfeReturn::Wolfe:
    return "Wolfe";
  case WolfeReturn::Armijo:
    return "Armijo";
  case WolfeReturn::ConvergedGradient:
    return "ConvergedGradient";
  case WolfeReturn::ConvergedObjective:
    return "ConvergedObjective";
  case WolfeReturn::ConvergedObjectiveAndGradient:
    return "ConvergedObjectiveAndGradient";
  case WolfeReturn::IntervalTooSmall:
    return "IntervalTooSmall";
  case WolfeReturn::StepTooSmall:
    return "StepTooSmall";
  case WolfeReturn::ReachedMaxStep:
    return "ReachedMaxStep";
  case WolfeReturn::NumericalIssue:
    return "NumericalIssue";
  case WolfeReturn::Fail:
    return "Fail";
  default:
    return "UNKNOWN";
  }
}

struct Eval {
  double alpha;   // alpha
  double obj;   // obj
  double dir;   // directional derivative
};

struct WolfeData {
  Eigen::VectorXd theta_;
  Eigen::VectorXd theta_grad_;
  Eigen::VectorXd a_;
  double obj_;
  double alpha_;
  double dir_;
  WolfeData(Eigen::Index n) :
    theta_(Eigen::VectorXd::Zero(n)),
    theta_grad_(Eigen::VectorXd::Zero(n)),
    a_(Eigen::VectorXd::Zero(n)),
    obj_(0.0),
    alpha_(0.0),
    dir_(0.0) {}
  void swap(WolfeData& other) {
    theta_.swap(other.theta_);
    theta_grad_.swap(other.theta_grad_);
    a_.swap(other.a_);
    std::swap(obj_, other.obj_);
    std::swap(alpha_, other.alpha_);
    std::swap(dir_, other.dir_);
  }
  void swap(WolfeData& other, Eval& eval) {
    theta_.swap(other.theta_);
    a_.swap(other.a_);
    std::swap(obj_, eval.obj);
    std::swap(alpha_, eval.alpha);
    std::swap(dir_, eval.dir);
  }
};
struct WolfeInfo {
  WolfeData curr_;
  WolfeData prev_;
  WolfeData scratch_;
  WolfeInfo(Eigen::Index n) :
    curr_(n),
    prev_(n),
    scratch_(n) {}
  void swap() {
    curr_.swap(prev_);
    curr_.swap(scratch_);
  }
};

/**
 * @brief  Strong‑Wolfe line search with cubic‑interpolation "zoom" for
 *Laplace‐style log‑likelihood problems.
 *
 * This routine searches along the space of the latent gaussian *a*
 * \f$a(\alpha) = prev.a_ + \alpha p`, `p = a − prev.a_\f$,
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
 *  void swap(WolfeData& other) {
    theta_.swap(other.theta_);
    theta_grad_.swap(other.theta_grad_);
    a_.swap(other.a_);
    std::swap(obj_, other.obj_);
    std::swap(alpha_, other.alpha_);
    std::swap(dir_, other.dir_);
  }

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
 *`opt.min_alpha`.
 *
 * * **Gradient reuse** – the caller provides `theta_grad` computed at the
 *   starting point; inside the loop it is overwritten with fresh gradients
 *   via `laplace_likelihood::theta_grad`.
 * * **Early‑exit for \f$\alpha\f$ = 1 or < `min_alpha`** Laplace problems
 *commonly accept a full Newton step. The function short‑circuits to avoid any
 *extra work.
 *
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
inline WolfeStatus wolfe_line_search(
    WolfeInfo& wolfe_info, F&& ll_fun, Obj&& obj_fun, Grad&& grad_fun,
    const Eigen::MatrixXd& covariance, LLArgs&& ll_args, Options&& opt,
    Stream* msgs) {
  auto& curr = wolfe_info.curr_;
  auto& prev = wolfe_info.prev_;
  auto& scratch = wolfe_info.scratch_;
  Eigen::VectorXd p = curr.a_ - prev.a_;
  double dir_deriv_init = grad_fun(prev.a_, prev.theta_, prev.theta_grad_).dot(p);
  Eval low{0.0, prev.obj_, dir_deriv_init};
  prev.dir_ = dir_deriv_init;
  auto armijo_ok = [&](const Eval& eval) -> bool {
    return check_armijo(eval.obj, prev.obj_, eval.alpha, dir_deriv_init, opt);
  };
  auto wolfe_ok = [&](const Eval& eval) -> bool {
    return check_wolfe_curve(eval.dir, dir_deriv_init, opt);
  };
  int total_updates = 0;
  auto update_step
      = [&p, &prev, &covariance, &ll_fun, &ll_args, &obj_fun, &grad_fun, msgs,
         &total_updates](auto& a_in, auto& theta_in, auto& theta_grad_in,
                         auto& eval_in) {
          total_updates++;
          a_in = prev.a_ + eval_in.alpha * p;
          theta_in = covariance * a_in;
          theta_grad_in
              = laplace_likelihood::theta_grad(ll_fun, theta_in, ll_args, msgs);
          eval_in.obj = obj_fun(a_in, theta_in);
          eval_in.dir = grad_fun(a_in, theta_in, theta_grad_in).dot(p);
        };
  double alpha_start = std::clamp(curr.alpha_ * 2.0, opt.min_alpha,
                                 opt.max_alpha);
  Eval high{alpha_start, curr.obj_, dir_deriv_init};
  update_step(scratch.a_, scratch.theta_, curr.theta_grad_, high);
  Eval best = low;  // keep the best Armijo-OK in case strong-Wolfe fails
  {
    while (!(std::isfinite(high.obj) && scratch.theta_.allFinite())) {
      high.alpha *= 0.5;
      if (high.alpha < opt.min_alpha) {
        debug::print("Exit on precheck numerical trouble", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::StepTooSmall, total_updates, 0};
      }
      update_step(scratch.a_, scratch.theta_, curr.theta_grad_, high);
    }
    debug::print("First precheck: ", 1, "high.alpha: ", high.alpha,
                 "high.obj:   ", high.obj, "deriv_high: ", high.dir,
                 "deriv_init: ", dir_deriv_init);
    // Quick accept if Armijo and Wolfe conditions are satisfied
    if (armijo_ok(high)) {
      if (wolfe_ok(high)) {
        debug::print("Exit on first precheck", 1);
        curr.swap(scratch, high);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, 0};
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
  //    "scratch.theta_:     ", scratch.theta_.transpose().eval());
  int loop_iter = 0;
  const auto grad_tol = opt.abs_grad_threshold;
  const auto obj_tol = opt.abs_obj_threshold;
  // If true we have already found a good first point
  bool found_right = false;
  int num_backtracks = 0;
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
  while (!found_right && high.alpha < opt.max_alpha) {
    num_backtracks++;
    // 1. Evaluate f(α) and g(α)
    update_step(scratch.a_, scratch.theta_, curr.theta_grad_, high);
    debug::print("First While", 1, "Second Iter:       ", loop_iter++,
                 "high.alpha: ", high.alpha, "high.obj:   ", high.obj,
                 "deriv_high: ", high.dir, "deriv_init: ", dir_deriv_init,
                 "scratch.theta_:  ", scratch.theta_.transpose());
    const bool finite_ok = std::isfinite(high.obj) && scratch.theta_.allFinite();
    // 2. Handle numerical trouble first
    if (!finite_ok) {  //   f or g is NaN/Inf → shrink
      high.alpha *= 0.5;
      if (high.alpha < opt.min_alpha) {
        break;
      }
      continue;
    }
    if (armijo_ok(high)) {
      // [1]
      if (wolfe_ok(high)) {
        curr.swap(scratch, high);
        debug::print("Exit on first while", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks};
      } else {
        if (best.obj < high.obj) {
          best = high;
        }
        // [2]
        if (std::signbit(high.dir)) {
          low = high;
          high.alpha *= opt.scale_up;
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
  auto check_bounds = [&](auto&& curr_eval) {
    // Check for grad convergence
    if (std::abs(curr_eval.dir) <= grad_tol ||  // tiny slope
      std::abs(curr_eval.obj - curr.obj_) <= obj_tol
          && curr_eval.alpha < 1e-8) {                // tiny gain
      if (curr_eval.obj != low.obj && armijo_ok(curr_eval)) {
        update_step(scratch.a_, scratch.theta_, curr.theta_grad_, curr_eval);
        curr.swap(scratch, curr_eval);
      }
      debug::print("total_updates", total_updates);
      if (std::abs(curr_eval.dir) <= grad_tol &&  // tiny slope
          std::abs(curr_eval.obj - curr.obj_) <= obj_tol) {
        debug::print("Exit on grad_tol and obj_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedObjectiveAndGradient, total_updates, num_backtracks};
      } else if (std::abs(curr_eval.dir) <= grad_tol) {
        debug::print("Exit on grad_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedGradient, total_updates, num_backtracks};
      } else if (std::abs(curr_eval.obj - curr.obj_) <= obj_tol) {
        debug::print("Exit on obj_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedObjective, total_updates, num_backtracks};
      } else {
        debug::print("Exit on alpha failure", 1);
        return WolfeStatus{WolfeReturn::IntervalTooSmall, total_updates, num_backtracks};
      }
    }
    return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks};
  };
  auto check_b = check_bounds(high);
  if (check_b.stop_ != WolfeReturn::Wolfe) {
    return check_b;
  }
  update_step(scratch.a_, scratch.theta_, curr.theta_grad_, high);

  debug::print("_______End First While: ", 1, "high.alpha: ", high.alpha,
               "high.obj:   ", high.obj, "high.dir: ", high.dir,
               "dir_deriv_init: ", dir_deriv_init);
  loop_iter = 0;
  Eval mid{low};
  while (high.alpha - low.alpha > opt.min_alpha) {
    num_backtracks++;
    const double diff_alpha = high.alpha - low.alpha;
    mid.alpha = cubic_or_bisect_max(low, high, opt);
    update_step(scratch.a_, scratch.theta_, curr.theta_grad_, mid);
    debug::print("Cube: ", 1, "Cube Iter:           ", loop_iter++,
                 "mid.alpha:      ", mid.alpha, "mid.obj:        ", mid.obj,
                 "high.dir: ", mid.dir,
                 "low.alpha:      ", low.alpha, "low.obj:        ", low.obj,
                 "low.dir:  ", low.dir,
                 "high.alpha:     ", high.alpha, "high.obj:       ", high.obj,
                 "high.dir: ", high.dir);
    const bool finite_ok = std::isfinite(mid.obj) && scratch.theta_.allFinite();
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
        curr.swap(scratch, mid);
        debug::print("Exit on safe on zoom", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks};
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
    auto check_bb = check_bounds(mid);
    if (check_bb.stop_ != WolfeReturn::Wolfe) {
      return check_bb;
    }

  }
  debug::print("Failed zoom: ", 1, "Failed zoom:", 1, "low.alpha: ", low.alpha,
               "low.obj:   ", low.obj, "deriv_low: ", low.dir,
               "high.alpha:", high.alpha, "high.obj:  ", high.obj,
               "deriv_high:", high.dir);
  // On failure, use the best point we have found so far that at least satisfies armijo
  update_step(scratch.a_, scratch.theta_, curr.theta_grad_, best);
  const bool armijo_ok_mid = armijo_ok(best);
  const bool curve_ok_mid = wolfe_ok(best);
  curr.a_.swap(scratch.a_);
  curr.theta_.swap(scratch.theta_);
  curr.theta_grad_ = laplace_likelihood::theta_grad(ll_fun, scratch.theta_, ll_args, msgs);
  // We already calculated best.obj and theta_grad so no need to recompute here
  best.dir = grad_fun(curr.a_, curr.theta_, curr.theta_grad_).dot(p);
  curr.obj_ = best.obj;
  curr.alpha_ = best.alpha;
  curr.dir_ = mid.dir;
  if (armijo_ok_mid && curve_ok_mid) {
    debug::print("Exit on safe after zoom", 1);
    debug::print("total_updates", total_updates);
    return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks};
  } else if (armijo_ok_mid) {
    debug::print("Exit on only satisfying armijo", 1);
    debug::print("total_updates", total_updates);
    return WolfeStatus{WolfeReturn::Armijo, total_updates, num_backtracks};
  } else {
    debug::print("Exit on failure", 1);
    debug::print("total_updates", total_updates);
    return WolfeStatus{WolfeReturn::Fail, total_updates, num_backtracks};
  }
}
}

}
#endif
