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
#include <type_traits>

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

/**
 * Selects a safeguarded trial point for maximizing a scalar function on a line.
 *
 * The routine assumes a 1-D bracket [x_left, x_right], together with function
 * values and directional derivatives at both endpoints. Internally it:
 *
 *   1. Normalizes the interval to $s \in [0, 1]$ via
 *        `x(s) = x_left + s * (x_right - x_left)`,
 *      and builds a cubic Hermite model F(s) that matches
 *      `{f_left, df_left}` at `s = 0` and `{f_right, df_right}` at `s = 1`.
 *
 *   2. Initializes the best candidate at the bisection point s = 0.5.
 *
 *   3. Adds a secant candidate for the derivative root, using the derivative
 *      values at s = 0 and s = 1:
 *        F'(0) = width * df_left, F'(1) = width * df_right.
 *
 *   4. Finds stationary points of the cubic model by solving F'(s) = 0.
 *      This reduces to a quadratic equation, handled with a numerically
 *      stable q-formula and tolerances for degeneracies (nearly linear
 *      derivative, small discriminant, etc.).
 *
 *   5. Evaluates all admissible model-based candidates and keeps the one
 *      with the largest F(s). All such candidates are restricted to a
 *      trimmed interior
 *        $s \in [edge_guard, 1 - edge_guard]$,
 *      i.e. $x \in [x_left + edge_guard * width, x_right - edge_guard *
 * width]$.
 *
 *   6. If the bracket is invalid (x_right <= x_left), any input is non-finite,
 *      or the interval is too tiny to be useful, it falls back to pure
 *      bisection and returns (x_left + x_right) / 2.
 *
 * The intended use is in line searches for MAXIMIZATION with a "well-formed"
 * bracket satisfying x_left < x_right, df_left > 0, df_right < 0. These sign
 * conditions are not required for safety; they only improve the model.
 *
 * @tparam Scalar  Floating-point scalar type (float, double, long double).
 *
 * @param x_left   Left endpoint of the current bracket.
 * @param f_left   Function value at x_left, i.e. f(x_left).
 * @param df_left  Directional derivative at x_left with respect to increasing
 * x, i.e. f'(x_left) in the search direction.
 * @param x_right  Right endpoint of the current bracket.
 * @param f_right  Function value at x_right, i.e. f(x_right).
 * @param df_right Directional derivative at x_right with respect to increasing
 * x, i.e. f'(x_right) in the search direction.
 *
 * @return A trial point in the trimmed interior of (x_left, x_right) chosen by
 *         the cubic/derivative model. If inputs are degenerate, the midpoint
 *         (x_left + x_right) / 2 is returned instead.
 */
template <typename Scalar>
[[nodiscard]] inline Scalar cubic_spline(Scalar x_left, Scalar f_left,
                                         Scalar df_left, Scalar x_right,
                                         Scalar f_right,
                                         Scalar df_right) noexcept {
  const Scalar midpoint = (x_left + x_right) / Scalar(2);

  // Basic validation: ordering + finiteness.
  if (!(x_right > x_left) || !std::isfinite(f_left) || !std::isfinite(f_right)
      || !std::isfinite(df_left) || !std::isfinite(df_right)) {
    return midpoint;
  }

  const Scalar width = x_right - x_left;
  const Scalar eps = std::numeric_limits<Scalar>::epsilon();

  // If the bracket is extremely tight, just bisect.
  {
    const Scalar x_scale
        = std::max(std::max(std::abs(x_left), std::abs(x_right)), Scalar(1));
    if (width <= eps * x_scale) {
      return midpoint;
    }
  }

  // Derivatives with respect to s, where x = x_left + s * width.
  const Scalar df_left_s = width * df_left;    // F'(0)
  const Scalar df_right_s = width * df_right;  // F'(1)

  // Cubic Hermite coefficients in $s \in [0,1]$:
  //   F(s) = a3*s^3 + a2*s^2 + a1*s + a0
  // with F(0) = f_left, F'(0) = df_left_s, F(1) = f_right, F'(1) = df_right_s.
  const Scalar a0 = f_left;
  const Scalar a1 = df_left_s;
  const Scalar a2
      = Scalar(3) * (f_right - f_left) - Scalar(2) * df_left_s - df_right_s;
  const Scalar a3 = Scalar(2) * (f_left - f_right) + df_left_s + df_right_s;

  auto eval = [&](Scalar s) -> Scalar {
    // Horner evaluation of F(s).
    return ((a3 * s + a2) * s + a1) * s + a0;
  };

  // Candidates are restricted to a trimmed interior [edge_guard, 1 -
  // edge_guard].
  constexpr Scalar edge_guard = Scalar(1e-9);

  Scalar best_s = 0.5;
  Scalar best_val = eval(0.5);
  auto consider = [&](Scalar s) {
    if (!std::isfinite(s)) {
      return;
    }
    if (!(s > edge_guard && s < Scalar(1) - edge_guard)) {
      return;
    }
    const Scalar value = eval(s);
    if (value > best_val) {
      best_s = s;
      best_val = value;
    }
  };

  // 1) Secant estimate for the derivative root between s = 0 and s = 1.
  {
    const Scalar denom = df_left_s - df_right_s;
    const Scalar deriv_scale = std::max(
        std::max(std::abs(df_left_s), std::abs(df_right_s)), Scalar(1));
    if (std::abs(denom) > eps * deriv_scale) {
      const Scalar s_secant
          = df_left_s / denom;  // Root of linear interpolation of F'.
      consider(s_secant);
    }
  }

  // 2) Stationary points of the cubic model (F'(s) = 0).
  {
    // F'(s) = 3*a3*s^2 + 2*a2*s + a1 = 0.
    const Scalar A = Scalar(3) * a3;
    const Scalar B = Scalar(2) * a2;
    const Scalar C = a1;

    const Scalar scale
        = std::max(std::max(std::abs(B), std::abs(C)), Scalar(1));
    const Scalar A_tol = eps * scale;

    if (std::abs(A) <= A_tol) {
      // Degenerate to (approximately) linear: B*s + C = 0.
      const Scalar B_tol = eps * scale;
      if (std::abs(B) > B_tol) {
        consider(-C / B);
      }
    } else {
      // Proper quadratic: A*s^2 + B*s + C = 0.
      Scalar disc = std::fma(-Scalar(4) * A, C, B * B);  // B^2 - 4AC
      const Scalar disc_scale
          = std::max(B * B + std::abs(Scalar(4) * A * C), Scalar(1));
      const Scalar disc_tol = Scalar(10) * eps * disc_scale;

      // Treat tiny negative discriminants as zero.
      if (disc < Scalar(0) && -disc <= disc_tol) {
        disc = Scalar(0);
      }

      if (disc >= Scalar(0)) {
        const Scalar r = std::sqrt(disc);
        const Scalar q = -Scalar(0.5) * (B + std::copysign(r, B));
        const Scalar q_scale = std::max(std::abs(B) + r, Scalar(1));
        const Scalar q_tol = eps * q_scale;

        if (std::abs(q) > q_tol) {
          const Scalar s1 = q / A;
          const Scalar s2 = C / q;
          consider(s1);
          consider(s2);
        } else {
          // Fallback: vertex of the quadratic derivative.
          const Scalar s_vertex = -B / (Scalar(2) * A);
          consider(s_vertex);
        }
      }
    }
  }

  return x_left + best_s * width;
}

template <typename Eval, typename Options>
inline auto cubic_spline(Eval&& low, Eval&& high, Options&& opt) {
  auto alpha = cubic_spline(low.alpha(), low.obj(), low.dir(), high.alpha(),
                            high.obj(), high.dir());
  const double width = high.alpha() - low.alpha();
  const double guard = 1e-3 * width;  // or make this an option
  alpha = std::clamp(alpha, low.alpha() + guard, high.alpha() - guard);
  return alpha;
}

template <typename Option>
inline auto check_armijo(double obj_next, double obj_init, double alpha_next,
                         double dir0, Option&& opt) {
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
  // When a check passes and we want to continue searching
  Continue
};

/**
 * Struct to hold the result status of the Wolfe line search.
 */
struct WolfeStatus {
  // total updates/evaluations
  int num_evals_{0};
  int num_backtracks_{-1};
  WolfeReturn stop_{WolfeReturn::Fail};
  // Whether a valid new step was found
  bool accept_{false};
  WolfeStatus() = default;
  WolfeStatus(WolfeReturn stop, int evals, int back)
      : num_evals_(evals), num_backtracks_(back), stop_(stop), accept_{false} {}
  WolfeStatus(WolfeReturn stop, int evals, int back, bool success)
      : num_evals_(evals),
        num_backtracks_(back),
        stop_(stop),
        accept_{success} {}
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
  inline auto& alpha() { return alpha_; }
  inline const auto& alpha() const { return alpha_; }
  inline auto& obj() { return obj_; }
  inline const auto& obj() const { return obj_; }
  inline auto& dir() { return dir_; }
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
  void update(WolfeData& other, const Eval& eval) {
    theta_.swap(other.theta_);
    a_.swap(other.a_);
    theta_grad_.swap(other.theta_grad_);
    eval_ = eval;
  }
  inline auto& theta() & { return theta_; }
  inline auto&& theta() && { return std::move(theta_); }
  inline const auto& theta() const& { return theta_; }

  inline auto& theta_grad() & { return theta_grad_; }
  inline const auto& theta_grad() const& { return theta_grad_; }
  inline auto&& theta_grad() && { return std::move(theta_grad_); }
  inline auto& a() & { return a_; }
  inline const auto& a() const& { return a_; }
  inline auto&& a() && { return std::move(a_); }
  inline auto& obj() & { return eval_.obj(); }
  inline const auto& obj() const& { return eval_.obj(); }
  inline auto& alpha() & { return eval_.alpha(); }
  inline const auto& alpha() const& { return eval_.alpha(); }
  inline auto& dir() & { return eval_.dir(); }
  inline const auto& dir() const& { return eval_.dir(); }
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
  /**
   * Construct WolfeInfo with a consistent (a_init, theta_init) pair.
   *
   * When the caller supplies a non-zero theta_init, the corresponding
   * a_init = Sigma^{-1} * theta_init must be provided so that the
   * invariant theta = Sigma * a holds at initialization.  This avoids
   * an inflated initial objective (the prior term -0.5 * a'*theta would
   * otherwise vanish when a is zero but theta is not).
   */
  template <typename ObjFun, typename Theta0, typename ThetaGradF>
  WolfeInfo(ObjFun&& obj_fun, const Eigen::VectorXd& a_init, Theta0&& theta0,
            ThetaGradF&& theta_grad_f, int /*tag*/)
      : curr_(std::forward<ObjFun>(obj_fun), a_init,
              std::forward<Theta0>(theta0),
              std::forward<ThetaGradF>(theta_grad_f)),
        prev_(curr_),
        scratch_(a_init.size()) {
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

  inline void flip_direction() {
    if (this->init_dir_ < 0) {
      this->p_ *= -1;
      this->init_dir_ *= -1;
    }
  }

  inline auto& curr() & { return curr_; }
  inline const auto& curr() const& { return curr_; }
  inline auto&& curr() && { return std::move(curr_); }
  inline auto& prev() & { return prev_; }
  inline const auto& prev() const& { return prev_; }
  inline auto&& prev() && { return std::move(prev_); }
  inline auto& scratch() & { return scratch_; }
  inline const auto& scratch() const& { return scratch_; }
  inline auto&& scratch() && { return std::move(scratch_); }
};

/**
 * Retry evaluation of a step until it passes a validity check.
 *
 * The update callable is invoked with `(curr, prev, eval, p)` and is expected
 * to fill `eval` (at `eval.alpha()`) with the objective and directional
 * derivative. If the evaluation is not valid, the backoff callable should
 * shrink `eval.alpha()` and return whether another retry should be attempted.
 * The validity check can inspect the evaluation and, for non-void updates, the
 * returned status.
 *
 * @tparam Update Callable that performs one evaluation step. Must accept 4
 * arguments.
 * @tparam Proposal Proposed step type passed to `update`.
 * @tparam Curr Current state type passed to `update`.
 * @tparam Prev Previous state type passed to `update`.
 * @tparam Eval Evaluation record containing alpha/obj/dir.
 * @tparam P Search direction type passed to `update`.
 * @tparam Backoff Callable that shrinks `eval.alpha()` and returns a bool.
 *
 * @param[in] update Evaluator invoked as `update(curr, prev, eval, p)`.
 * @param[in,out] proposal Proposed step forwarded to `update`.
 * @param[in] curr Current state forwarded to `update`.
 * @param[in] prev Previous state forwarded to `update`.
 * @param[in,out] eval Evaluation record, updated in-place by `update`.
 * @param[in] p Search direction forwarded to `update`.
 * @param[in] backoff Shrinks alpha and returns whether another retry should
 * occur.
 *
 * @return For void updates, returns void. Otherwise returns the value from the
 *         first valid evaluation.
 */
template <typename Update, typename Proposal, typename Curr, typename Prev,
          typename Eval, typename P, typename Backoff>
inline auto retry_evaluate(Update&& update, Proposal&& proposal, Curr&& curr,
                           Prev&& prev, Eval& eval, P&& p, Backoff&& backoff) {
  while (true) {
    auto res = update(proposal, curr, prev, eval, p);
    if (res) {
      return res;
    }
    if (!backoff(eval)) {
      return res;
    }
  }
}

/**
 * @brief Strong Wolfe line search for maximization.
 *
 * Finds step $\alpha$ along $a(\alpha) = a_0 + \alpha\,p$ satisfying the
 * strong Wolfe conditions for **maximization**:
 *
 * - Armijo (sufficient increase):
 *   $\phi(\alpha) \ge \phi(0) + c_1\,\alpha\,\phi'(0)$
 * - Curvature:
 *   $|\phi'(\alpha)| \le c_2\,|\phi'(0)|$
 *
 * where $\phi(\alpha)$ is the objective and
 * $\phi'(\alpha) = \nabla\phi^\top p$ is the directional derivative.
 * Falls back to Armijo-only or convergence statuses when strong Wolfe
 * cannot be achieved.
 *
 * The search maintains a left endpoint `low`, a right endpoint `high`,
 * and a fallback `best` (best Armijo-satisfying point seen so far).
 * Non-finite evaluations are contracted by factor `opt.tau` automatically.
 *
 * ## Phase 1: Aggresive Expansion
 * The user given initial step size is expanded by `opt.scale_up` until
 * either Wolfe conditions are violated or `opt.max_alpha` is reached.
 * If the first evaluation satisfies both conditions, the search expands
 * further ("zoom-up") to find the largest such step before accepting.
 *
 * ## Phase 2: Expansion / Bracketing
 *
 * Starting from $\alpha_0 = \text{clamp}(\text{curr.alpha} \cdot
 * \text{scale\_up},\; \text{min\_alpha},\; \text{max\_alpha})$,
 * evaluates increasing step sizes to find a bracket `[low, high]`
 * guaranteed to contain a Wolfe point. If the first evaluation already
 * satisfies both conditions, the search expands further ("zoom-up")
 * to find the largest such step before accepting.
 *
 * ```
 * | Armijo | Wolfe | f'    | Action                    |
 * |--------|-------|-------|---------------------------|
 * | T      | T     | -     | Accept (return Wolfe)     |
 * | T      | F     | > 0   | low = high, expand high   |
 * | T      | F     | <= 0  | Bracket found, go to zoom |
 * | F      | -     | -     | Bracket found, go to zoom |
 * ```
 *
 * The goal of this phase is to find a valid bracket `[low, high]`
 * Ideally such that the low endpoint satisfies Armijo and has
 * a positive derivative, while the high endpoint has a negative
 * derivative. This gives us a shape like /\ to search through
 * in the zoom phase. If such a bracket cannot be found
 * within the step size limits, that is fine but the zoom phase
 * will start with bisections instead of the cubic interpolation.
 *
 * ## Phase 3: Zoom
 *
 * Narrows `[low, high]` to find a strong Wolfe point. Uses cubic
 * interpolation when the bracket has a derivative sign change and
 * `high` satisfies Armijo; otherwise bisects.
 *
 * ```
 * | Armijo | f(mid)>f(lo) | Wolfe | f'(mid) | Action     |
 * |--------|--------------|-------|---------|------------|
 * | T      | T            | T     | -       | Accept     |
 * | T      | T            | F     | > 0     | low = mid  |
 * | T      | T            | F     | <= 0    | high = mid |
 * | T      | F            | -     | -       | high = mid |
 * | F      | -            | -     | -       | high = mid |
 * ```
 *
 * If no strong Wolfe point is found, falls back to the best
 * Armijo-satisfying point (status `Armijo`) or `Fail` if none exists.
 *
 * ## UpdateFun contract
 *
 * ```
 * void update_fun(WolfeData& out, WolfeData& curr,
 *                 WolfeData& prev, Eval& e, const P& p);
 * ```
 *
 * Must compute $\phi(\alpha)$ and $\phi'(\alpha)$ at `e.alpha()`,
 * store results in `e.obj()` and `e.dir()`, and populate `out` with
 * the full candidate state (theta, theta_grad, a). Must not mutate
 * `curr` or `prev`.
 *
 * @tparam Info      Wolfe working set (prev_, curr_, scratch_, p_, init_dir_)
 * @tparam UpdateFun Step evaluator satisfying the contract above
 * @tparam Options   Provides c1, c2, tau, scale_up, min/max_alpha,
 *                   max_iterations, and convergence thresholds
 * @tparam Stream    Output stream type (may be nullptr)
 *
 * @param wolfe_info Working set with baseline (prev_), current (curr_),
 *                   scratch space, direction p_, and initial derivative
 * @param update_fun Step evaluator callback
 * @param opt        Line search options
 * @param msgs       Optional diagnostic stream (nullptr to disable)
 *
 * @return `WolfeStatus` with:
 *  - `Wolfe` -- both conditions satisfied, `curr_` updated
 *  - `Armijo` -- only Armijo satisfied, `curr_` updated
 *  - `ConvergedGradient` / `ConvergedObjective` /
 *    `ConvergedObjectiveAndGradient` -- early convergence,
 *    `curr_` updated iff Armijo holds
 *  - `IntervalTooSmall` -- bracket collapsed, `curr_` updated iff Armijo holds
 *  - `StepTooSmall` -- $\alpha$ contracted below `min_alpha`
 *  - `ReachedMaxStep` -- evaluation budget exceeded
 *  - `Fail` -- no Armijo point found, `curr_` unchanged
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
  int total_updates = 0;
  Eval best = low;  // keep the best Armijo-OK in case strong-Wolfe fails
  auto update_with_tick = [&total_updates, &opt, &best, &update_fun](
                              auto&& proposal, auto&& curr, auto&& prev,
                              Eval& e, auto&& p) {
    const bool over_budget = total_updates > opt.max_iterations;
    if (over_budget) {
      // Soft budget: stop evaluating new trial points once exceeded.
      if (check_armijo(best, prev, opt)) {
        update_fun(proposal, curr, prev, best, p);
        curr.update(proposal, best);
        if (check_wolfe(best, prev, opt)) {
          return WolfeStatus{WolfeReturn::Wolfe, total_updates, 0, true};
        } else {
          return WolfeStatus{WolfeReturn::Armijo, total_updates, 0, true};
        }
      }
      return WolfeStatus{WolfeReturn::ReachedMaxStep, total_updates, 0, false};
    } else {
      update_fun(proposal, curr, prev, e, p);
      ++total_updates;
      return WolfeStatus{WolfeReturn::Continue, total_updates, 0, false};
    }
  };
  double alpha_start
      = std::clamp(curr.alpha() * opt.scale_up, opt.min_alpha, opt.max_alpha);
  Eval high{alpha_start, curr.obj(), dir_deriv_init};
  WolfeStatus wolfe_check{WolfeReturn::Continue, 0, 0, false};
  // Initial check for numerical trouble
  {
    wolfe_check = update_with_tick(scratch, curr, prev, high, p);
    if (wolfe_check.stop_ != WolfeReturn::Continue) {
      return wolfe_check;
    }
    if (high.alpha() < opt.min_alpha) {
      return WolfeStatus{WolfeReturn::StepTooSmall, total_updates, 0, false};
    }
    // Quick accept if Armijo and Wolfe conditions are satisfied
    if (check_armijo(high, prev, opt)) {
      if (check_wolfe(high, prev, opt)) {
        // Try zooming up till we hit a fail
        best = high;
        while (check_armijo(high, prev, opt) && check_wolfe(high, prev, opt)) {
          best = high;
          high.alpha() *= opt.scale_up;
          if (high.alpha() > opt.max_alpha) {
            break;
          }
          wolfe_check = update_with_tick(scratch, curr, prev, high, p);
          if (wolfe_check.stop_ != WolfeReturn::Continue) {
            return wolfe_check;
          }
        }
        wolfe_check = update_with_tick(scratch, curr, prev, best, p);
        if (wolfe_check.stop_ != WolfeReturn::Continue) {
          return wolfe_check;
        }
        curr.update(scratch, best);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, 0, true};
      } else {
        if (best.obj() < high.obj()) {
          best = high;
        }
      }
    }
  }
  int num_backtracks = 0;
  /**
   * For each case
   * | armijo     | wolfe | sign(g) | Action
   * -------+-------+---------+--------------------------------
   * | [1]  T     |   T   |         | Accept alpha
   * | [2]  T     |   F   |   > 0   | set low=high, expand high
   * | [3]  T     |   F   |   < 0   | stop
   * | [4]  F     |   T   |         | stop
   * | [5]  F     |   F   |         | stop
   **/
  while (high.alpha() < opt.max_alpha) {
    num_backtracks++;
    // 1. Evaluate f(alpha) and g(alpha)
    wolfe_check = update_with_tick(scratch, curr, prev, high, p);
    if (wolfe_check.stop_ != WolfeReturn::Continue) {
      return wolfe_check;
    }
    if (check_armijo(high, prev, opt)) {
      if (check_wolfe(high, prev, opt)) {  // [1]
        curr.update(scratch, high);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks,
                           true};
      }
      if (best.obj() < high.obj()) {
        best = high;
      }
      if (high.dir() > 0) {  // [2]
        low = high;
        high.alpha() *= opt.scale_up;
        continue;
      }
    }
    // [3,4,5]
    break;
  }
  const double grad_tol
      = std::max(opt.abs_grad_threshold,
                 opt.rel_grad_threshold * std::abs(dir_deriv_init));
  const double obj_tol
      = std::max(opt.abs_obj_threshold,
                 opt.rel_obj_threshold * (1.0 + std::abs(prev.obj())));
  auto check_bounds = [&](auto&& curr_eval) {
    // Check for grad convergence
    const bool slope_check = std::abs(curr_eval.dir()) <= grad_tol;
    // tiny slope or gain
    const bool obj_check = std::abs(curr_eval.obj() - prev.obj()) <= obj_tol;
    // alpha too small
    const bool alpha_check = curr_eval.alpha() < opt.min_alpha;
    if (slope_check || obj_check || alpha_check) {
      bool step_ok
          = curr_eval.obj() != low.obj() && check_armijo(curr_eval, prev, opt);
      if (slope_check && obj_check) {
        return WolfeStatus{WolfeReturn::ConvergedObjectiveAndGradient,
                           total_updates, num_backtracks, step_ok};
      } else if (slope_check) {
        return WolfeStatus{WolfeReturn::ConvergedGradient, total_updates,
                           num_backtracks, step_ok};
      } else if (obj_check) {
        return WolfeStatus{WolfeReturn::ConvergedObjective, total_updates,
                           num_backtracks, step_ok};
      } else {
        return WolfeStatus{WolfeReturn::IntervalTooSmall, total_updates,
                           num_backtracks, step_ok};
      }
    }
    return WolfeStatus{WolfeReturn::Continue, total_updates, num_backtracks,
                       false};
  };
  auto check_b = check_bounds(high);
  if (check_b.stop_ != WolfeReturn::Continue) {
    if (check_b.accept_) {
      curr.update(scratch, high);
    }
    return check_b;
  }
  wolfe_check = update_with_tick(scratch, curr, prev, high, p);
  if (wolfe_check.stop_ != WolfeReturn::Continue) {
    return wolfe_check;
  }
  /**
   *  Zoom phase
   * | #  | Armijo | f(m) > f(l) | Wolfe | f'(mid) |  Action    |
   * |----|--------|-------------|-------|---------|------------|
   * | 1  |   T    |      T      |   T   |    -    | Accept     |
   * | 2  |   T    |      T      |   F   |   > 0   | low = mid  |
   * | 3  |   T    |      T      |   F   |   <= 0  | high = mid |
   * | 4  |   T    |      F      |   -   |    -    | high = mid |
   * | 5  |   F    |      -      |   -   |    -    | high = mid |
   */
  while ((high.alpha() - low.alpha() > opt.min_alpha)
         && high.alpha() > opt.min_alpha) {
    num_backtracks++;
    const bool have_sign_change = (low.dir() * high.dir() < 0);
    const bool high_armijo_ok = check_armijo(high, prev, opt);
    const bool use_cubic = have_sign_change && high_armijo_ok;
    // Choose trial alpha: cubic when bracket is good, else bisection.
    double alpha_mid{0};
    if (use_cubic) {
      alpha_mid = cubic_spline(low, high, opt);
    } else {
      alpha_mid = 0.5 * (low.alpha() + high.alpha());
    }
    if (alpha_mid <= opt.min_alpha) {
      break;
    }
    Eval mid{alpha_mid, 0.0, 0.0};
    auto wolfe_check = update_with_tick(scratch, curr, prev, mid, p);
    if (wolfe_check.stop_ != WolfeReturn::Continue) {
      return wolfe_check;
    }
    if (mid.alpha() <= opt.min_alpha) {
      return WolfeStatus{WolfeReturn::StepTooSmall, total_updates,
                         num_backtracks, false};
    }
    if (check_armijo(mid, prev, opt)) {
      if (check_wolfe(mid, prev, opt)) {  // [1]
        curr.update(scratch, mid);
        return WolfeStatus{WolfeReturn::Wolfe, total_updates, num_backtracks,
                           true};
      }
      // Track best Armijo-OK point for fallback.
      if (mid.obj() > best.obj()) {
        best = mid;
      }
      if (mid.obj() > low.obj()) {
        if (mid.dir() > 0) {  // [2]
          low = mid;
        } else {  // [3]
          high = mid;
        }
      } else {
        // [4]
        high = mid;
      }
    } else {
      // [5]
      high = mid;
    }
    // Convergence/guard-rail checks (uses prev/grad_tol/obj_tol etc.)
    auto bounds_check = check_bounds(mid);
    if (bounds_check.stop_ != WolfeReturn::Continue) {
      if (bounds_check.accept_) {
        curr.update(scratch, mid);
      }
      return bounds_check;
    }
  }
  // On failure, use the best point we have found so far that at least satisfies
  // armijo
  const bool armijo_ok_best = check_armijo(best, prev, opt);
  if (armijo_ok_best) {
    wolfe_check = update_with_tick(scratch, curr, prev, best, p);
    curr.update(scratch, best);
    return WolfeStatus{WolfeReturn::Armijo, total_updates, num_backtracks,
                       true};
  } else {
    return WolfeStatus{WolfeReturn::Fail, total_updates, num_backtracks, false};
  }
}
}  // namespace internal

}  // namespace stan::math
#endif
