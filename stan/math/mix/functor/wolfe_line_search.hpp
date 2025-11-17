#ifndef STAN_MATH_MIX_FUNCTOR_WOLFE_LINE_SEARCH_HPP
#define STAN_MATH_MIX_FUNCTOR_WOLFE_LINE_SEARCH_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
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
  constexpr explicit laplace_line_search_options(int max_iter)
      : max_iterations(max_iter) {}
  constexpr laplace_line_search_options() = default;
  int max_iterations{1000};
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

  double rel_grad_threshold{1e-3};  // off by default
  double rel_obj_threshold{1e-10};  // off by default
};
namespace internal {
namespace debug {

constexpr void print(int tabs) {}

#ifdef LAPLACE_DEBUG
std::ofstream csv_file("../wolfe_line_search_debug.csv");
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
inline Scalar cubic_or_bisect_max(Scalar a, Scalar fa, Scalar fpa, Scalar b,
                                  Scalar fb, Scalar fpb) {
  // Basic validation.
  if (!(b > a) || !std::isfinite(fa) || !std::isfinite(fb)
      || !std::isfinite(fpa) || !std::isfinite(fpb)) {
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
  struct Candidate {
    Scalar s;
    Scalar val;
  };
  Candidate best{Scalar(0.5), eval(Scalar(0.5))};  // start with bisection

  auto consider = [&](Scalar s) {
    if (!(s > Scalar(0) && s < Scalar(1)) || !std::isfinite(s))
      return;
    const Scalar v = eval(s);
    if (v > best.val)
      best = {s, v};
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
        const Scalar s2 = (q == Scalar(0))
                              ? std::numeric_limits<Scalar>::quiet_NaN()
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
  constexpr Scalar hi = Scalar(0.0) * (1.0 - (Scalar(1) - edge_guard))
                        + (1.0 * (Scalar(1) - edge_guard));
  const Scalar s_star = std::clamp(best.s, lo, hi);
  // Map back to alpha-space.
  return a + s_star * width;
}

template <typename Eval, typename Options>
inline auto cubic_or_bisect_max(Eval&& low, Eval&& high, Options&& opt) {
  auto alpha = cubic_or_bisect_max(low.alpha(), low.obj(), low.dir(),
                                   high.alpha(), high.obj(), high.dir());
  const double width = high.alpha() - low.alpha();
  const double guard = 1e-3 * width;  // or make this an option
  alpha = std::clamp(alpha, low.alpha() + guard, high.alpha() - guard);
  return alpha;
}

template <typename Option>
inline auto check_armijo(double obj_next, double obj_init, double alpha_next,
                         double dir0, Option&& opt) {
  debug::print(
      "check_armijo: ", 2, "armijo:    ",
      (obj_next >= obj_init + alpha_next * dir0 * opt.c1 ? 1 : 0),
      "obj_next:   ", obj_next, "obj_init:   ", obj_init,
      "alpha_next: ", alpha_next, "dir0:       ", dir0, "c1:         ", opt.c1,
      "obj + alpha * dir0 * c1: ", (obj_init + alpha_next * dir0 * opt.c1));
  return (obj_next >= obj_init)
         && (obj_next >= obj_init + alpha_next * dir0 * opt.c1);
}

template <typename Eval, typename WolfeT, typename Option>
inline bool check_armijo(const Eval& eval, const WolfeT& prev,
                         const Option& opt) {
  return check_armijo(eval.obj(), prev.obj(), eval.alpha(), prev.dir(), opt);
}

template <typename Option>
inline auto check_wolfe_curve(double dir_deriv_next, double dir_deriv_init,
                              Option&& opt) {
  debug::print(
      "check_wolfe_curve: ", 2, "wolfe:    ",
      (std::abs(dir_deriv_next) <= (opt.c2 * std::abs(dir_deriv_init)) ? 1 : 0),
      "deriv_next: ", dir_deriv_next, "deriv_init: ", dir_deriv_init,
      "c2:         ", opt.c2, "abs(d_next):   ", std::abs(dir_deriv_next),
      "abs(d_init)*c2 ", (std::abs(dir_deriv_init) * opt.c2));
  return std::abs(dir_deriv_next) <= (opt.c2 * std::abs(dir_deriv_init));
}

template <typename Eval, typename WolfeT, typename Option>
inline bool check_wolfe(const Eval& eval, const WolfeT& prev,
                        const Option& opt) {
  return check_wolfe_curve(eval.dir(), prev.dir(), opt);
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
  Fail,
  Continue
};

/**
 * Struct to hold the result status of the Wolfe line search.
 */
struct WolfeStatus {
  // total updates/evaluations
  int num_evals_{0};
  int num_backtracks_{0};
  WolfeReturn stop_{WolfeReturn::Fail};
  // Whether a valid new step was found
  bool success_{false};
  WolfeStatus() = default;
  WolfeStatus(WolfeReturn stop, int evals, int back)
      : num_evals_(evals),
        num_backtracks_(back),
        stop_(stop),
        success_{false} {}
  WolfeStatus(WolfeReturn stop, int evals, int back, bool success)
      : num_evals_(evals),
        num_backtracks_(back),
        stop_(stop),
        success_{success} {}
};

/**
 * Helper function for pretty printing
 */
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
    case WolfeReturn::Continue:
      return "Continue";
    default:
      return "UNKNOWN";
  }
}

/**
 * evaluation struct for Wolfe line search
 */
struct Eval {
  // alpha
  double alpha_{0.0};
  // obj
  double obj_{0.0};
  // directional derivative
  double dir_{0.0};
  inline auto&& alpha() { return alpha_; }
  inline const auto& alpha() const { return alpha_; }
  inline auto&& obj() { return obj_; }
  inline const auto& obj() const { return obj_; }
  inline auto&& dir() { return dir_; }
  inline const auto& dir() const { return dir_; }
  constexpr Eval(double alpha, double obj, double dir)
      : alpha_(alpha), obj_(obj), dir_(dir) {}
  constexpr Eval() = default;
};

/**
 * Data used in current evaluation of wolfe line search at a particular stepsize
 */
struct WolfeData {
  // current parameter values
  Eigen::VectorXd theta_;
  // parameter gradients
  Eigen::VectorXd theta_grad_;
  // latent variable
  Eigen::VectorXd a_;
  // evaluation data
  Eval eval_;
  explicit WolfeData(Eigen::Index n)
      : theta_(Eigen::VectorXd::Zero(n)),
        theta_grad_(Eigen::VectorXd::Zero(n)),
        a_(Eigen::VectorXd::Zero(n)),
        eval_() {}

  template <typename ObjFun, typename ThetaGradF, typename Theta0>
  WolfeData(ObjFun&& obj_fun, const Eigen::VectorXd& a, const Theta0& theta0,
            ThetaGradF&& theta_grad_f)
      : theta_(theta0),
        theta_grad_(theta_grad_f(theta_)),
        a_(a),
        eval_(1.0, obj_fun(a_, theta_), 0.0) {}

  template <typename ObjFun, typename ThetaGradF, typename Theta0>
  WolfeData(ObjFun&& obj_fun, Eigen::Index n, const Theta0& theta0,
            ThetaGradF&& theta_grad_f)
      : theta_(theta0),
        theta_grad_(theta_grad_f(theta_)),
        a_(Eigen::VectorXd::Zero(n)),
        eval_(1.0, obj_fun(a_, theta_), 0.0) {}

  // Initialize with theta = 0
  template <typename LLFun, typename LLArgs, typename Msgs>
  WolfeData(Eigen::Index n, const LLFun& ll_fun, const LLArgs& ll_args,
            const Msgs& msgs)
      : WolfeData(n, Eigen::VectorXd::Zero(n), ll_fun, ll_args, msgs) {}

  void update(WolfeData& other) {
    theta_.swap(other.theta_);
    theta_grad_.swap(other.theta_grad_);
    a_.swap(other.a_);
    eval_ = other.eval_;
  }
  void update(WolfeData& other, Eval& eval) {
    theta_.swap(other.theta_);
    a_.swap(other.a_);
    theta_grad_.swap(other.theta_grad_);
    eval_ = eval;
  }
  inline auto&& theta() { return theta_; }
  inline const auto& theta() const { return theta_; }
  inline auto&& theta_grad() { return theta_grad_; }
  inline const auto& theta_grad() const { return theta_grad_; }
  inline auto&& a() { return a_; }
  inline const auto& a() const { return a_; }
  inline auto&& obj() { return eval_.obj(); }
  inline const auto& obj() const { return eval_.obj(); }
  inline auto&& alpha() { return eval_.alpha(); }
  inline const auto& alpha() const { return eval_.alpha(); }
  inline auto&& dir() { return eval_.dir(); }
  inline const auto& dir() const { return eval_.dir(); }
};

/**
 * Data object used in wolfe line search
 */
struct WolfeInfo {
  // Current step data. On output will be the proposal step
  WolfeData curr_;
  // Previous step data
  WolfeData prev_;
  // Scratch space for evaluations of proposal steps
  WolfeData scratch_;
  // Search direction
  Eigen::VectorXd p_;
  // Initial directional derivative
  double init_dir_;
  template <typename ObjFun, typename Theta0, typename ThetaGradF>
  WolfeInfo(ObjFun&& obj_fun, Eigen::Index n, Theta0&& theta0,
            ThetaGradF&& theta_grad_f)
      : curr_(std::forward<ObjFun>(obj_fun), n, std::forward<Theta0>(theta0),
              std::forward<ThetaGradF>(theta_grad_f)),
        prev_(curr_),
        scratch_(n) {
    if (!std::isfinite(curr_.obj())) {
      throw std::domain_error(
          "laplace_marginal_density: log likelihood is not finite at initial "
          "theta and likelihood arguments.");
    }
  }
  WolfeInfo(WolfeData&& curr, WolfeData&& prev)
      : curr_(std::move(curr)),
        prev_(std::move(prev)),
        scratch_(curr_.theta().size()),
        p_(Eigen::VectorXd::Zero(curr_.theta().size())),
        init_dir_(0.0) {}
  explicit WolfeInfo(Eigen::Index n)
      : curr_(n),
        prev_(curr_),
        scratch_(n),
        p_(Eigen::VectorXd::Zero(n)),
        init_dir_(0.0) {}
};

/**
 * @brief Strong-Wolfe line search with expansion, bracketing, and
 * cubic/bisection zoom
 *
 * This routine searches along the ray
 * \f[
 *   a(\alpha) = a_0 + \alpha\,p,\qquad p = a_1 - a_0,
 * \f]
 * to find the largest step \f$\alpha\f$ that satisfies the **strong-Wolfe**
 * conditions
 *
 * \f{align*}{
 *   \phi(\alpha) &\le \phi(0) + c_1\,\alpha\,\phi'(0) \quad\text{(Armijo)},\\
 *   |\phi'(\alpha)| &\le c_2\,|\phi'(0)| \quad\text{(curvature)},
 * \f}
 *
 * where \f$\phi(\alpha)\f$ is the objective at \f$a(\alpha)\f$ and
 * \f$\phi'(\alpha)=\nabla\phi(a(\alpha))^\top p\f$ is the directional
 * derivative.
 *
 * ## How the search proceeds (phases)
 *
 * The algorithm works with a *left* point `low` and a *right* point `high`,
 * stores the “best so far” Armijo-OK point `best`, and reuses a scratch buffer
 * to avoid recomputation when a step is finally accepted.
 *
 * 1. **Numerical contraction (pre-check).**
 *    Start from \f$\alpha_0=\mathrm{clamp}( \text{curr.alpha} \cdot
 * \text{opt.scale_up}, \text{opt.min_alpha}, \text{opt.max_alpha})\f$. Evaluate
 * \f$\phi(\alpha_0)\f$ and the derived state. If either the objective or
 * derived quantities (e.g. \f$\theta(\alpha)\f$, its gradient) are non-finite,
 *    shrink \f$\alpha \leftarrow \alpha \cdot \text{opt.tau}\f$ until
 * everything is finite.
 *    **Ends when** values become finite → continue, or
 * \f$\alpha<\text{opt.min_alpha}\f$ → **returns** `StepTooSmall`.
 *
 * 2. **Quick accept & “zoom-up”.**
 *    If the current `high` satisfies **both** Armijo and curvature, repeatedly
 *    *expand* \f$\alpha \leftarrow \alpha \cdot \text{opt.scale_up}\f$ while
 * the strong-Wolfe tests continue to pass, keeping the last valid step as
 * `best`.
 *    **Ends when** a test fails, \f$\alpha>\text{opt.max_alpha}\f$, or
 * iteration limit is reached. The function **returns** the last passing step
 * with `Wolfe`.
 *
 * 3. **Expansion & right-bracketing.**
 *    Otherwise, continue evaluating larger \f$\alpha\f$ values:
 *    - If Armijo holds and \f$\phi'(\alpha)>0\f$, promote the left endpoint
 *      `low ← high` and *expand* again.
 *      (We are on the “far side” of the minimum.)
 *    - If Armijo holds but \f$\phi'(\alpha)\le 0\f$, or if Armijo fails, a
 * right endpoint has been found -> stop expanding and proceed to zoom.
 *    **Ends when** a right endpoint is found,
 * \f$\alpha>\text{opt.max_alpha}\f$, or iteration limit is reached.
 *
 * 4. **Safety bisection (optional).**
 *    If both endpoints currently have positive directional derivatives, take
 * pure bisection steps until a sign change is observed or the interval is
 * small.
 *
 * 5. **Cubic/bisection zoom.**
 *    Repeatedly propose an interior point \f$\alpha_m\f$ using a cubic
 * interpolation of the bracket end-points (falling back to bisection when the
 * cubic is degenerate), evaluate there, and shrink the bracket according to:
 *    - If Armijo and curvature hold → **returns** `Wolfe`.
 *    - If Armijo holds but curvature fails → keep the better side; prefer the
 * sub-interval with a derivative sign change.
 *    - If Armijo fails → move the right endpoint left.
 *    Throughout zoom we remember the best Armijo-OK point `best` in case
 *    strong-Wolfe cannot be satisfied.
 *    **Ends when** the interval width is \f$\le\text{opt.min_alpha}\f$, a
 * convergence tolerance triggers, or the iteration limit is reached.
 *
 * 6. **Convergence/guard rails.**
 *    After each evaluation we check:
 *    - Gradient tolerance: \f$|\phi'(\alpha)| \le
 *      \max(\text{opt.abs\_grad\_threshold},\
 * \text{opt.rel\_grad\_threshold}\,|\phi'(0)|)\f$.
 *    - Objective change tolerance (for very small steps):
 *      \f$|\phi(\alpha)-\phi(0)| \le
 *      \max(\text{opt.abs\_obj\_threshold},\
 * \text{opt.rel\_obj\_threshold}\,(1+|\phi(0)|))\f$ and \f$\alpha <
 * \text{opt.min_alpha}\f$. If triggered, the routine returns one of
 * `ConvergedGradient`, `ConvergedObjective`, or
 * `ConvergedObjectiveAndGradient`. If neither test is met but the bracket
 * becomes too small, it returns `IntervalTooSmall` (accepting the current point
 * only if it satisfies Armijo).
 *
 * If zoom exits without a strong-Wolfe point but an Armijo-OK point was seen,
 * the routine **returns** `Armijo` at that best point; otherwise it **returns**
 * `Fail`.
 *
 * ### Contract for the step evaluator `UpdateFun`
 *
 * The line search delegates *all* numerical work to a user-supplied callable
 * `update_fun` that evaluates the objective and its directional derivative at a
 * trial step and prepares the candidate state for acceptance.
 *
 * **Required signature**
 * @code
 * void update_fun(WolfeData& out, Eval& e, const Direction& p);
 * @endcode
 *
 * **Inputs**
 *  - `e.alpha()` is the trial step \f$\alpha\f$ to evaluate.
 *  - `p` is the search direction.
 *  - The baseline (point at \f$\alpha=0\f$) and its values live in
 *    `wolfe_info.prev_` (made available when you constructed/captured
 * `update_fun`).
 *
 * **You must:**
 *  - Form the trial point \f$a(\alpha)=a_0+\alpha p\f$ and any derived
 * quantities your model needs (e.g.,
 *  \f$\theta(\alpha)\f$, \f$\nabla\theta(\alpha)\f$).
 *  - Compute the objective \f$\phi(\alpha)\f$ and set `e.obj()`.
 *  - Compute the directional derivative
 * \f$\phi'(\alpha)=\nabla\phi(a(\alpha))^\top p\f$ and set `e.dir()`.
 *  - Populate `out` with the *entire* candidate state corresponding to
 * \f$\alpha\f$ (e.g., parameters, `theta()`, `theta_grad()`, cached factors,
 * etc.). The search will call `out.swap(...)` to accept this state; **do not**
 * mutate `curr_`/`prev_`.
 *
 * **Postconditions**
 *  - `e.obj()` and `e.dir()` are finite when the step is numerically valid.
 *  - `out.theta()` and `out.theta_grad()` are finite (used by the search to
 *    detect numerical trouble).
 *
 * `update_fun` is free to close over `ll_fun` and `ll_args`. This routine
 * itself does not call `ll_fun` directly; it simply invokes `update_fun`
 * whenever it needs an evaluation.
 *
 * ### Options and tolerances
 * The `opt` object must provide:
 *  - `double c1, c2;` // Wolfe constants (0 < c1 < c2 < 1)
 *  - `double tau;` // contraction factor (0 < tau < 1)
 *  - `double scale_up;` // expansion factor (> 1)
 *  - `double min_alpha, max_alpha;` // step bounds
 *  - `int    max_iterations;` // max number of calls to
 * `update_fun`
 *  - `double abs_grad_threshold, rel_grad_threshold;`
 *  - `double abs_obj_threshold,  rel_obj_threshold;`
 *
 * @tparam Info Type providing the fields documented under *Wolfe line search
 * data*.
 * @tparam UpdateFun Callable used to evaluate a trial step; must satisfy the
 * contract above.
 * @tparam Options Type providing the fields documented under *Options and
 * tolerances*.
 * @tparam Stream Output stream type for diagnostics (e.g., `std::ostream`); may
 * be `nullptr`.
 *
 * @param wolfe_info Line-search working set (holds `prev_`, `curr_`,
 * `scratch_`, the direction `p_`, and the initial directional derivative).
 * @param update_fun Callable that evaluates \f$\phi(\alpha)\f$,
 * \f$\phi'(\alpha)\f$, and prepares the candidate state (see contract above).
 * @param opt Wolfe and algorithmic options.
 * @param msgs Optional stream for debug messages; pass `nullptr` to disable.
 *
 * @return `WolfeStatus` describing how the search ended:
 *  - `Wolfe` — Found a step satisfying **both** Armijo and curvature; `curr_`
 * updated.
 *  - `Armijo` — Could not satisfy curvature; returns the best Armijo-OK step;
 * `curr_` updated.
 *  - `ConvergedGradient` — \f$|\phi'(\alpha)|\f$ below tolerance at a tiny
 * step; `curr_` may be updated iff that point also satisfies Armijo.
 *  - `ConvergedObjective` — Objective improvement below tolerance at a tiny
 * step; same acceptance rule.
 *  - `ConvergedObjectiveAndGradient` — Both tests met; same acceptance rule.
 *  - `IntervalTooSmall` — Bracket is smaller than `min_alpha` without meeting
 * strong-Wolfe; accepts the current point only if Armijo holds.
 *  - `StepTooSmall` — Numerical contraction drove \f$\alpha<\text{min_alpha}\f$
 * before a finite evaluation.
 *  - `ReachedMaxStep` — Exceeded `max_iterations`; if the best Armijo-OK point
 * exists it is returned (and may be `Wolfe` or `Armijo`), otherwise search
 * aborts without updating `curr_`.
 *  - `Fail` — No Armijo-OK point was found; `curr_` is left unchanged.
 *
 * The status also reports:
 *  - `num_evals` — total calls to `update_fun`.
 *  - `num_backtracks` — times the step was reduced during bracketing/zoom.
 *  - `accepted` — whether `curr_` was updated to a new point.
 */
template <typename Info, typename UpdateFun, typename Options, typename Stream>
inline WolfeStatus wolfe_line_search(Info& wolfe_info, UpdateFun&& update_fun,
                                     Options&& opt, Stream* msgs) {
  auto& curr = wolfe_info.curr_;
  auto& prev = wolfe_info.prev_;
  auto& scratch = wolfe_info.scratch_;
  auto&& p = wolfe_info.p_;
  auto&& dir_deriv_init = wolfe_info.init_dir_;
  Eval low{0.0, prev.obj(), dir_deriv_init};
  prev.dir() = dir_deriv_init;
  auto armijo_ok = [&prev, &opt](const Eval& eval) -> bool {
    return check_armijo(eval, prev, opt);
  };
  auto wolfe_ok = [&prev, &opt](const Eval& eval) -> bool {
    return check_wolfe(eval, prev, opt);
  };
  int total_updates = 0;
  auto assign_step
      = [](WolfeData& out, WolfeData& buf, Eval& e) { out.update(buf, e); };
  auto update_with_tick = [&total_updates, &update_fun, &prev, &curr](
                              WolfeData& buf, Eval& e, auto&& p) {
    update_fun(buf, curr, prev, e, p);
    ++total_updates;
  };
  auto check_max_steps = [&assign_step, &p, &total_updates, &armijo_ok,
                          &wolfe_ok, &update_with_tick,
                          &msgs](auto&& scratch, auto&& curr, auto&& prev,
                                 auto&& best, auto&& opt) {
    if (total_updates > opt.max_iterations) {
      debug::print("Exit on precheck max iterations", 1);
      debug::print("total_updates", total_updates);
      if (armijo_ok(best)) {
        update_with_tick(scratch, best, p);
        assign_step(curr, scratch, best);
        if (wolfe_ok(best)) {
          return WolfeStatus{WolfeReturn::Wolfe, total_updates, 0, true};
        } else {
          return WolfeStatus{WolfeReturn::Armijo, total_updates, 0, true};
        }
      }
      return WolfeStatus{WolfeReturn::ReachedMaxStep, total_updates, 0, false};
    } else {
      return WolfeStatus{WolfeReturn::Continue, total_updates, 0, false};
    }
  };
  double alpha_start
      = std::clamp(curr.alpha() * opt.scale_up, opt.min_alpha, opt.max_alpha);
  Eval high{alpha_start, curr.obj(), dir_deriv_init};
  update_with_tick(scratch, high, p);
  Eval best = low;  // keep the best Armijo-OK in case strong-Wolfe fails
  {
    while (!(std::isfinite(high.obj()) && scratch.theta().allFinite())) {
      high.alpha() *= opt.tau;
      if (high.alpha() < opt.min_alpha) {
        debug::print("Exit on precheck numerical trouble", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::StepTooSmall, total_updates, 0, false};
      }
      update_with_tick(scratch, high, p);
      auto check_steps = check_max_steps(scratch, curr, prev, high, opt);
      if (check_steps.stop_ != WolfeReturn::Continue) {
        return check_steps;
      }
    }
    debug::print("First precheck: ", 1, "high.alpha(): ", high.alpha(),
                 "high.obj():   ", high.obj(), "deriv_high: ", high.dir(),
                 "deriv_init: ", dir_deriv_init);
    // Quick accept if Armijo and Wolfe conditions are satisfied
    if (armijo_ok(high)) {
      if (wolfe_ok(high)) {
        // Try zooming up till we hit a fail
        best = high;
        while (armijo_ok(high) && wolfe_ok(high)) {
          best = high;
          auto check_steps = check_max_steps(scratch, curr, prev, best, opt);
          if (check_steps.stop_ != WolfeReturn::Continue) {
            return check_steps;
          }
          high.alpha() *= opt.scale_up;
          if (high.alpha() > opt.max_alpha) {
            break;
          }
          update_with_tick(scratch, high, p);
          debug::print("Zoom up: ", 1, "high.alpha(): ", high.alpha(),
                       "high.obj():   ", high.obj(), "deriv_high: ", high.dir(),
                       "deriv_init: ", dir_deriv_init);
        }
        debug::print("Exit on first precheck", 1);
        update_with_tick(scratch, best, p);
        assign_step(curr, scratch, best);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, 0, true};
      } else {
        if (best.obj() < high.obj()) {
          best = high;
        }
      }
    }
  }
  // If current alpha fails, backtrack down till we find a good point
  debug::print("Begin Loop: ", 1, "Initial alpha: ", high.alpha());
  int loop_iter = 0;
  const double grad_tol
      = std::max(opt.abs_grad_threshold,
                 opt.rel_grad_threshold * std::abs(dir_deriv_init));
  const double obj_tol
      = std::max(opt.abs_obj_threshold,
                 opt.rel_obj_threshold * (1.0 + std::abs(prev.obj())));
  // If true we have already found a good first point
  bool found_right = false;
  int num_backtracks = 0;
  /**
   * For each case
   * | armijo     | wolfe | sign(g) | Action
   * -------+-------+---------+--------------------------------
   * | [1]  T     |   T   |         | Accept alpha
   * | [2]  T     |   F   |   > 0   | set low=high, expand high
   * | [3]  T     |   F   |   < 0   | Set alpha_high <- alpha, stop
   * | [4]  F     |   T   |         | Set alpha_high <- alpha, stop
   * | [5]  F     |   F   |         | Set alpha_high <- alpha, stop
   **/
  while (!found_right && high.alpha() < opt.max_alpha) {
    num_backtracks++;
    // 1. Evaluate f(alpha) and g(alpha)
    update_with_tick(scratch, high, p);
    debug::print("First While", 1, "Second Iter:       ", loop_iter++,
                 "high.alpha(): ", high.alpha(), "high.obj():   ", high.obj(),
                 "deriv_high: ", high.dir(), "deriv_init: ", dir_deriv_init,
                 "scratch.theta():  ", scratch.theta().transpose());
    const bool finite_ok
        = std::isfinite(high.obj()) && scratch.theta().allFinite();
    // 2. Handle numerical trouble first
    if (!finite_ok) {  //   f or g is NaN/Inf → shrink
      high.alpha() *= 0.5;
      if (high.alpha() < opt.min_alpha) {
        break;
      }
      continue;
    }
    if (armijo_ok(high)) {
      // [1]
      if (wolfe_ok(high)) {
        assign_step(curr, scratch, high);
        debug::print("Exit on first while", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks,
                           true};
      } else {
        if (best.obj() < high.obj()) {
          best = high;
        }
        // [2]
        if (high.dir() > 0) {
          low = high;
          high.alpha() *= opt.scale_up;
          continue;
        } else {
          // [3]
          found_right = true;
        }
      }
    }
    // [4,5]
    found_right = true;
    auto check_steps = check_max_steps(scratch, curr, prev, best, opt);
    if (check_steps.stop_ != WolfeReturn::Continue) {
      return check_steps;
    }
  }
  auto check_bounds = [&](auto&& curr_eval) {
    // Check for grad convergence
    if (std::abs(curr_eval.dir()) <= grad_tol ||  // tiny slope
        std::abs(curr_eval.obj() - prev.obj()) <= obj_tol
            && curr_eval.alpha() < opt.min_alpha) {  // tiny gain
      bool step_ok = curr_eval.obj() != low.obj() && armijo_ok(curr_eval);
      if (step_ok) {
        update_with_tick(scratch, curr_eval, p);
        assign_step(curr, scratch, curr_eval);
      }
      debug::print("total_updates", total_updates);
      if (std::abs(curr_eval.dir()) <= grad_tol &&  // tiny slope
          std::abs(curr_eval.obj() - prev.obj()) <= obj_tol) {
        debug::print("Exit on grad_tol and obj_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedObjectiveAndGradient,
                           total_updates, num_backtracks, step_ok};
      } else if (std::abs(curr_eval.dir()) <= grad_tol) {
        debug::print("Exit on grad_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedGradient, total_updates,
                           num_backtracks, step_ok};
      } else if (std::abs(curr_eval.obj() - prev.obj()) <= obj_tol) {
        debug::print("Exit on obj_tol", 1);
        return WolfeStatus{WolfeReturn::ConvergedObjective, total_updates,
                           num_backtracks, step_ok};
      } else {
        debug::print("Exit on alpha failure", 1);
        return WolfeStatus{WolfeReturn::IntervalTooSmall, total_updates,
                           num_backtracks, step_ok};
      }
    }
    return WolfeStatus{WolfeReturn::Continue, total_updates, num_backtracks,
                       false};
  };
  auto check_b = check_bounds(high);
  if (check_b.stop_ != WolfeReturn::Continue) {
    return check_b;
  }
  update_with_tick(scratch, high, p);
  loop_iter = 0;
  // Pure bisection to try to get a sign change.
  while ((low.dir() > 0 && high.dir() > 0)
         && (high.alpha() - low.alpha() > opt.min_alpha)) {
    Eval mid{0.5 * (low.alpha() + high.alpha()), 0.0, 0.0};
    update_with_tick(scratch, mid, p);
    if (!std::isfinite(mid.obj()) || !std::isfinite(mid.dir())
        || !scratch.theta().allFinite() || !scratch.theta_grad().allFinite()) {
      high.alpha() *= opt.tau;
      if (high.alpha() < opt.min_alpha)
        break;
      continue;
    }
    if (best.obj() < mid.obj()) {
      best = mid;
    }
    if (mid.dir() > 0) {
      low = mid;
    } else {
      high = mid;
    }
  }
  Eval mid{high};
  while (mid.alpha() > opt.min_alpha) {
    num_backtracks++;
    const double diff_alpha = high.alpha() - low.alpha();
    mid.alpha() = cubic_or_bisect_max(low, high, opt);
    if (mid.alpha() <= opt.min_alpha) {
      break;
    }
    update_with_tick(scratch, mid, p);
    if (mid.alpha() <= opt.min_alpha) {
      break;
    }
    const bool finite_ok
        = std::isfinite(mid.obj()) && scratch.theta().allFinite();
    if (!finite_ok) {
      high = mid;
      continue;
    }
    /**
     * |Armijo|Wolfe|mid.obj > low.obj|mid.dir < 0|Action|
     * |  T   |  T  |                 |           | Accept |
     * |  T   |  F  |       T         |     T     | high = mid |
     * |  T   |  F  |       T         |     F     | low = mid |
     * |  T   |  F  |       F         |     F     | low = mid |
     * |  T   |  F  |       F         |     F     | high = mid |
     * |  F   |  F  |       F         |     F     | high = mid |
     */
    if (armijo_ok(mid)) {
      if (wolfe_ok(mid)) {
        assign_step(curr, scratch, mid);
        debug::print("Exit on safe on zoom", 1);
        debug::print("total_updates", total_updates);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks,
                           true};
      } else {
        if (best.obj() < mid.obj()) {
          best = mid;
        }
        if (mid.obj() > low.obj()) {
          // sign change
          if (mid.dir() * low.dir() < 0) {
            high = mid;
          } else {
            low = mid;
          }
        } else {
          high = mid;
        }
      }
    } else {
      if (best.obj() < high.obj()) {
        best = high;
      }
      high = mid;
    }
    auto check_bb = check_bounds(mid);
    if (check_bb.stop_ != WolfeReturn::Continue) {
      return check_bb;
    }
    auto check_steps = check_max_steps(scratch, curr, prev, best, opt);
    if (check_steps.stop_ != WolfeReturn::Continue) {
      return check_steps;
    }
  }
  debug::print("Failed zoom: ", 1, "Failed zoom:", 1,
               "low.alpha(): ", low.alpha(), "low.obj():   ", low.obj(),
               "deriv_low: ", low.dir(), "high.alpha():", high.alpha(),
               "high.obj():  ", high.obj(), "deriv_high:", high.dir());
  // On failure, use the best point we have found so far that at least satisfies
  // armijo
  const bool armijo_ok_best = armijo_ok(best);
  const bool curve_ok_best = wolfe_ok(best);
  if (armijo_ok_best) {
    update_with_tick(scratch, best, p);
    assign_step(curr, scratch, best);
    debug::print("Exit on only satisfying armijo", 1);
    debug::print("total_updates", total_updates);
    return WolfeStatus{WolfeReturn::Armijo, total_updates, num_backtracks,
                       true};
  } else {
    debug::print("Exit on failure", 1);
    debug::print("total_updates", total_updates);
    return WolfeStatus{WolfeReturn::Fail, total_updates, num_backtracks, false};
  }
}
}  // namespace internal

}  // namespace stan::math
#endif
