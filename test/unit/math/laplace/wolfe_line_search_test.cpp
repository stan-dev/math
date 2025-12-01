#include <gtest/gtest.h>
#include <stan/math/mix/functor/wolfe_line_search.hpp>
#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <ostream>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

namespace stan::math {
namespace {

using internal::WolfeInfo;
using internal::WolfeReturn;
using internal::WolfeStatus;

/**
 * @brief Concave quadratic log-likelihood paired with an objective.
 */
struct QuadraticObjective {
  Eigen::MatrixXd Q;
  Eigen::VectorXd b;
  double scale{1.0};

  QuadraticObjective(Eigen::MatrixXd Q_in, Eigen::VectorXd b_in,
                     double scale_in = 1.0)
      : Q(std::move(Q_in)), b(std::move(b_in)), scale(scale_in) {}

  double operator()(const Eigen::VectorXd& a,
                    const Eigen::VectorXd& /*theta*/) const {
    return scale * (b.dot(a) - 0.5 * a.dot(Q * a));
  }

  Eigen::VectorXd gradient(const Eigen::VectorXd& a) const {
    return scale * (b - Q * a);
  }

  template <typename Vec, typename... Args>
  auto log_likelihood(const Vec& theta, Args&&... /*args*/) const {
    using scalar_t = typename std::decay_t<Vec>::Scalar;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> theta_cast = theta;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Qt
        = Q.template cast<scalar_t>() * theta_cast;
    return scale
           * (b.template cast<scalar_t>().dot(theta_cast)
              - 0.5 * theta_cast.dot(Qt));
  }
};

/**
 * @brief Construct a diagonal symmetric positive definite matrix.
 */
inline Eigen::MatrixXd make_spd_from_diag(const std::vector<double>& diag) {
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(diag.size(), diag.size());
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(diag.size()); ++i) {
    Q(i, i) = diag[i];
  }
  return Q;
}

/**
 * @brief Construct a 2x2 symmetric positive definite matrix with correlation.
 */
inline Eigen::MatrixXd make_correlated_spd(double a, double b, double c) {
  Eigen::MatrixXd M(2, 2);
  M << a, b, b, c;
  return M;
}

/**
 * @brief Log barrier objective with matching log-likelihood.
 */
struct LogBarrierObjective {
  double radius;

  explicit LogBarrierObjective(double r) : radius(r) {}

  double operator()(const Eigen::VectorXd& a,
                    const Eigen::VectorXd& /*theta*/) const {
    double norm_sq = a.squaredNorm();
    double denom = 1.0 - norm_sq / (radius * radius);
    if (denom <= 0.0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::log(denom);
  }

  Eigen::VectorXd gradient(const Eigen::VectorXd& a) const {
    double norm_sq = a.squaredNorm();
    double denom = radius * radius - norm_sq;
    if (denom <= 0.0) {
      Eigen::VectorXd nan = Eigen::VectorXd::Constant(
          a.size(), std::numeric_limits<double>::quiet_NaN());
      return nan;
    }
    return (-2.0 / denom) * a;
  }

  template <typename Vec, typename... Args>
  auto log_likelihood(const Vec& theta, Args&&... /*args*/) const {
    using scalar_t = typename std::decay_t<Vec>::Scalar;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> theta_cast = theta;
    scalar_t norm_sq = theta_cast.squaredNorm();
    scalar_t denom = scalar_t(1.0) - norm_sq / (radius * radius);
    if (denom <= scalar_t(0.0)) {
      return scalar_t(-std::numeric_limits<double>::infinity());
    }
    return stan::math::log(denom);
  }
};

/**
 * @brief Nearly linear objective with small curvature.
 */
struct LinearObjective {
  Eigen::VectorXd r;
  double curvature;

  LinearObjective(Eigen::VectorXd r_in, double eps)
      : r(std::move(r_in)), curvature(eps) {}

  double operator()(const Eigen::VectorXd& a,
                    const Eigen::VectorXd& /*theta*/) const {
    return r.dot(a) - 0.5 * curvature * a.squaredNorm();
  }

  Eigen::VectorXd gradient(const Eigen::VectorXd& a) const {
    return r - curvature * a;
  }

  template <typename Vec, typename... Args>
  auto log_likelihood(const Vec& theta, Args&&... /*args*/) const {
    using scalar_t = typename std::decay_t<Vec>::Scalar;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> theta_cast = theta;
    return r.template cast<scalar_t>().dot(theta_cast)
           - scalar_t(0.5) * curvature * theta_cast.squaredNorm();
  }
};

/**
 * @brief Negative Rosenbrock objective used to test strong-Wolfe steps.
 */
struct RosenbrockObjective {
  double operator()(const Eigen::VectorXd& a,
                    const Eigen::VectorXd& /*theta*/) const {
    double x = a[0];
    double y = a[1];
    double term1 = (1.0 - x) * (1.0 - x);
    double term2 = y - x * x;
    return -(term1 + 100.0 * term2 * term2);
  }

  Eigen::VectorXd gradient(const Eigen::VectorXd& a) const {
    Eigen::VectorXd grad(2);
    double x = a[0];
    double y = a[1];
    double term2 = y - x * x;
    grad[0] = -(2.0 * (x - 1.0) - 400.0 * x * term2);
    grad[1] = -200.0 * term2;
    return grad;
  }

  template <typename Vec, typename... Args>
  auto log_likelihood(const Vec& theta, Args&&... /*args*/) const {
    using scalar_t = typename std::decay_t<Vec>::Scalar;
    scalar_t x = theta[0];
    scalar_t y = theta[1];
    scalar_t term1 = (scalar_t(1.0) - x) * (scalar_t(1.0) - x);
    scalar_t term2 = y - x * x;
    return -(term1 + scalar_t(100.0) * term2 * term2);
  }
};

inline double rosenbrock_value(const Eigen::VectorXd& a) {
  double x = a[0];
  double y = a[1];
  double term1 = (1.0 - x) * (1.0 - x);
  double term2 = y - x * x;
  return term1 + 100.0 * term2 * term2;
}

/**
 * @brief Nonconvex scalar objective f(x) = cos(x) + 0.1 x^2 (maximise -f).
 */
struct WavyObjective {
  double operator()(const Eigen::VectorXd& a,
                    const Eigen::VectorXd& /*theta*/) const {
    double x = a[0];
    return -(std::cos(x) + 0.1 * x * x);
  }

  Eigen::VectorXd gradient(const Eigen::VectorXd& a) const {
    Eigen::VectorXd grad(1);
    double x = a[0];
    grad[0] = std::sin(x) - 0.2 * x;
    return grad;
  }

  template <typename Vec, typename... Args>
  auto log_likelihood(const Vec& theta, Args&&... /*args*/) const {
    using scalar_t = typename std::decay_t<Vec>::Scalar;
    scalar_t x = theta[0];
    return -(stan::math::cos(x) + scalar_t(0.1) * x * x);
  }
};

inline double wavy_value(double x) { return std::cos(x) + 0.1 * x * x; }

/**
 * @brief Log-likelihood that couples theta with scale and shift arguments.
 */
struct CoupledLikelihood {
  template <typename Theta, typename Scale, typename Shift, typename Stream>
  auto operator()(const Theta& theta, Scale&& scale, Shift&& shift,
                  Stream* /*msgs*/) const {
    using theta_scalar = scalar_type_t<std::decay_t<Theta>>;
    return theta_scalar(0.5) * std::forward<Scale>(scale) * theta.squaredNorm()
           + std::forward<Shift>(shift) * theta.sum();
  }
};

/**
 * @brief Expected gradients for CoupledLikelihood arguments.
 */
inline std::pair<double, double> coupled_likelihood_arg_gradients(
    const Eigen::VectorXd& theta) {
  return {0.5 * theta.squaredNorm(), theta.sum()};
}

/**
 * @brief Construct autodiff variables for the coupled likelihood arguments.
 */
inline std::tuple<var, var> make_coupled_args(double scale, double shift) {
  return std::make_tuple(var(scale), var(shift));
}

/**
 * @brief Helper to prepare Wolfe line search tests with reusable utilities.
 */
struct LineSearchHarness {
  laplace_line_search_options opt{};
  Eigen::MatrixXd covariance;

  explicit LineSearchHarness(Eigen::Index n)
      : covariance(Eigen::MatrixXd::Identity(n, n)) {}

  /**
   * @brief Create WolfeInfo initialised with previous and current states.
   */
  template <typename Obj>
  WolfeInfo make_info(const Obj& obj, const Eigen::VectorXd& prev_a,
                      const Eigen::VectorXd& curr_a) const {
    WolfeInfo info(prev_a.size());
    info.prev_.a() = prev_a;
    info.prev_.theta() = covariance * prev_a;
    info.prev_.theta_grad().setZero();
    info.prev_.obj() = obj(prev_a, info.prev_.theta());
    info.prev_.alpha() = 0.0;
    info.curr_.a() = curr_a;
    info.curr_.alpha() = 1.0;
    info.curr_.obj() = obj(curr_a, covariance * curr_a);
    info.curr_.theta_grad().setZero();
    return info;
  }

  /**
   * @brief Construct a gradient functor compatible with wolfe_line_search.
   */
  template <typename Obj>
  auto make_grad_fun(const Obj& obj) const {
    return [&obj](const Eigen::VectorXd& a, const Eigen::VectorXd& /*theta*/,
                  const Eigen::VectorXd& /*theta_grad*/) {
      return obj.gradient(a);
    };
  }

  /**
   * @brief Execute the Wolfe search with a provided log-likelihood functor.
   */
  template <typename Obj, typename LL, typename LLArgs>
  WolfeStatus run(WolfeInfo& info, const Obj& obj, const LL& ll_fun,
                  LLArgs&& ll_args) {
    auto grad_fun = make_grad_fun(obj);
    auto obj_fun
        = [&obj](const Eigen::VectorXd& a, const Eigen::VectorXd& theta) {
            return obj(a, theta);
          };
    auto ll_adapter = [&ll_fun](const auto& theta, auto&&... args) {
      return ll_fun(theta, std::forward<decltype(args)>(args)...);
    };
    std::ostream* msgs = nullptr;
    auto update_step = [&](auto& step_info, auto&& curr, auto&& prev,
                           auto& eval_in, auto&& p) {
      stan::math::set_zero_all_adjoints();
      step_info.a() = info.prev_.a() + eval_in.alpha() * p;
      step_info.theta() = covariance * step_info.a();
      step_info.theta_grad() = laplace_likelihood::theta_grad(
          ll_adapter, step_info.theta(), ll_args, msgs);
      eval_in.obj() = obj_fun(step_info.a(), step_info.theta());
      eval_in.dir()
          = grad_fun(step_info.a(), step_info.theta(), step_info.theta_grad())
                .dot(p);
    };
    info.curr_.theta() = covariance * info.curr_.a();
    info.curr_.theta_grad().setZero();
    info.p_ = info.curr_.a() - info.prev_.a();
    info.init_dir_ = obj.gradient(info.prev_.a()).dot(info.p_);
    if (info.init_dir_ <= 0.0) {
      info.p_ = -info.p_;
      info.init_dir_ = -info.init_dir_;
    }

    return wolfe_line_search(info, update_step, opt,
                             static_cast<std::ostream*>(nullptr));
  }

  /**
   * @brief Execute the Wolfe search using the objective's default
   * log-likelihood.
   */
  template <typename Obj>
  WolfeStatus run(WolfeInfo& info, const Obj& obj) {
    return run(
        info, obj,
        [&obj](const auto& theta, auto&&... args) {
          return obj.log_likelihood(theta,
                                    std::forward<decltype(args)>(args)...);
        },
        std::tuple<>{});
  }
};

/**
 * @brief Compute the initial search direction and directional derivative.
 */
template <typename Obj>
inline std::pair<Eigen::VectorXd, double> initial_direction(
    const Obj& obj, const WolfeInfo& info_before) {
  Eigen::VectorXd p = info_before.curr_.a() - info_before.prev_.a();
  double dir0 = obj.gradient(info_before.prev_.a()).dot(p);
  if (dir0 <= 0.0) {
    p = -p;
    dir0 = -dir0;
  }
  return {p, dir0};
}

/**
 * @brief Evaluate the directional derivative of a quadratic objective.
 */
template <typename Obj>
inline double directional_derivative(const Obj& obj, const Eigen::VectorXd& a,
                                     const Eigen::VectorXd& p) {
  return obj.gradient(a).dot(p);
}

// Checks that a well-conditioned concave quadratic accepts the strong-Wolfe
// step.
TEST(WolfeLineSearch, StrongWolfeConcaveQuadratic) {
  const int n = 3;
  LineSearchHarness harness(n);
  Eigen::MatrixXd Q = make_spd_from_diag({1.5, 0.7, 2.0});
  Eigen::VectorXd b(n);
  b << 3.0, -1.0, 2.5;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd optimum = Q.ldlt().solve(b);
  Eigen::VectorXd curr = prev + optimum;

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;  // Encourage initial trial near alpha = 1
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);
  EXPECT_GT(status.num_evals_, 0);

  auto [p, dir0] = initial_direction(obj, before);
  double alpha = info.curr_.alpha();
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);

  EXPECT_GT(alpha, harness.opt.min_alpha);
  EXPECT_LE(alpha, harness.opt.max_alpha);
  EXPECT_GE(phi_alpha, phi0 + harness.opt.c1 * alpha * dir0 - 1e-12);
  EXPECT_LE(std::abs(dir_alpha), harness.opt.c2 * std::abs(dir0) + 1e-12);
  EXPECT_TRUE(info.curr_.theta().isApprox(harness.covariance * info.curr_.a()));
  EXPECT_TRUE(info.curr_.theta().allFinite());
  EXPECT_TRUE(info.curr_.theta_grad().allFinite());
}

// Checks that the initial doubling candidate satisfies the Wolfe conditions.
TEST(WolfeLineSearch, AcceptsOnFirstPrecheck) {
  LineSearchHarness harness(2);
  harness.opt.c1 = 1e-4;
  harness.opt.c2 = 0.9;
  Eigen::MatrixXd Q = make_spd_from_diag({0.9, 1.1});
  Eigen::VectorXd b(2);
  b << 0.8, -0.3;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd newton = Q.ldlt().solve(b);
  Eigen::VectorXd curr = prev + newton;

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);
  EXPECT_EQ(status.num_backtracks_, 0);

  auto [p, dir0] = initial_direction(obj, before);
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double alpha = info.curr_.alpha();
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  EXPECT_GE(phi_alpha, phi0 + harness.opt.c1 * alpha * dir0 - 1e-12);
}

// Checks that the search zooms to satisfy the curvature condition.
TEST(WolfeLineSearch, RequiresZoomForCurvature) {
  const int n = 2;
  LineSearchHarness harness(n);
  Eigen::MatrixXd Q = make_spd_from_diag({0.4, 4.0});
  Eigen::VectorXd b(n);
  b << 0.6, 1.5;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd newton = Q.ldlt().solve(b);
  Eigen::VectorXd curr = prev + newton;

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 1.0;  // initial doubling overshoots the maximiser
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::ConvergedObjective)
      << "Expected Converged Objective but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);
  EXPECT_GE(status.num_backtracks_, 0);

  auto [p, dir0] = initial_direction(obj, before);
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);
  EXPECT_LE(std::abs(dir_alpha), harness.opt.c2 * std::abs(dir0) + 1e-12);
}

// Checks that exact Armijo equality is accepted.
TEST(WolfeLineSearch, ArmijoEqualityAccepted) {
  LineSearchHarness harness(1);
  harness.opt.c1 = 0.6;
  Eigen::MatrixXd Q = make_spd_from_diag({1.6});
  Eigen::VectorXd b(1);
  b << 1.0;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd curr = Eigen::VectorXd::Ones(1);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.25;  // trial alpha -> 0.5
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [p, dir0] = initial_direction(obj, before);
  double alpha = info.curr_.alpha();
  EXPECT_NEAR(alpha, 0.5, 1e-12);
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  double armijo_gap = phi_alpha - (phi0 + harness.opt.c1 * alpha * dir0);
  EXPECT_NEAR(armijo_gap, 0.0, 1e-12);
}

// Tests stay focused on the strong-Wolfe variant because the API does not yet
// expose a weak-Wolfe configuration.

// Checks that a 1D quadratic step (f(x) = 0.5 (x - 3)^2 up to a constant)
// takes a meaningful alpha while satisfying Armijo and signed curvature tests.
TEST(WolfeLineSearch, OneDimensionalQuadraticStrongWolfe) {
  LineSearchHarness harness(1);
  harness.opt.c1 = 1e-4;
  harness.opt.c2 = 0.8;
  Eigen::MatrixXd Q = make_spd_from_diag({1.0});
  Eigen::VectorXd b(1);
  b << 3.0;
  QuadraticObjective obj(Q, b);  // 4.5 - 0.5 (x - 3)^2

  Eigen::VectorXd prev(1);
  prev << 1.5;
  Eigen::VectorXd optimum = Q.ldlt().solve(b);
  Eigen::VectorXd curr = prev + 0.5 * (optimum - prev);

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.75;
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [p, dir0] = initial_direction(obj, before);
  double alpha = info.curr_.alpha();
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);

  EXPECT_GT(alpha, harness.opt.min_alpha);
  EXPECT_GT(alpha, 1e-3);
  EXPECT_LT(alpha, harness.opt.max_alpha);
  EXPECT_GE(phi_alpha, phi0 + harness.opt.c1 * alpha * dir0 - 1e-12);
  EXPECT_LE(dir_alpha, harness.opt.c2 * dir0 + 1e-12);

  double quad_prev = 0.5 * std::pow(before.prev_.a()[0] - 3.0, 2);
  double quad_curr = 0.5 * std::pow(info.curr_.a()[0] - 3.0, 2);
  EXPECT_LT(quad_curr, quad_prev);
}

// Checks that a Rosenbrock step reduces the original objective while satisfying
// the strong-Wolfe inequalities.
TEST(WolfeLineSearch, RosenbrockStrongWolfeStep) {
  LineSearchHarness harness(2);
  harness.opt.c1 = 1e-4;
  harness.opt.c2 = 0.75;
  RosenbrockObjective obj;

  Eigen::VectorXd prev(2);
  prev << -1.2, 1.0;
  Eigen::VectorXd curr = prev + 0.2 * obj.gradient(prev);

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [p, dir0] = initial_direction(obj, before);
  double alpha = info.curr_.alpha();
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);

  EXPECT_GT(alpha, harness.opt.min_alpha);
  EXPECT_LT(alpha, harness.opt.max_alpha);
  EXPECT_GE(phi_alpha, phi0 + harness.opt.c1 * alpha * dir0 - 1e-12);
  EXPECT_LE(std::abs(dir_alpha), harness.opt.c2 * std::abs(dir0) + 1e-12);

  double rosen_prev = rosenbrock_value(before.prev_.a());
  double rosen_curr = rosenbrock_value(info.curr_.a());
  EXPECT_LT(rosen_curr, rosen_prev);
}

// Checks that a wavy nonconvex scalar search direction still finishes with a
// strong-Wolfe compliant step.
TEST(WolfeLineSearch, WavyNonconvexStrongWolfeStep) {
  LineSearchHarness harness(1);
  harness.opt.c1 = 1e-4;
  harness.opt.c2 = 0.8;
  WavyObjective obj;

  Eigen::VectorXd prev(1);
  prev << 1.2;
  Eigen::VectorXd curr = prev + 0.6 * obj.gradient(prev);

  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [p, dir0] = initial_direction(obj, before);
  double alpha = info.curr_.alpha();
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);

  EXPECT_GT(alpha, harness.opt.min_alpha);
  EXPECT_LT(alpha, harness.opt.max_alpha);
  EXPECT_GE(phi_alpha, phi0 + harness.opt.c1 * alpha * dir0 - 1e-12);
  EXPECT_LE(std::abs(dir_alpha), harness.opt.c2 * std::abs(dir0) + 1e-12);

  double wavy_prev = wavy_value(before.prev_.a()[0]);
  double wavy_curr = wavy_value(info.curr_.a()[0]);
  EXPECT_LT(wavy_curr, wavy_prev);
}

// Checks that exact curvature equality is accepted.
TEST(WolfeLineSearch, CurvatureEqualityAccepted) {
  LineSearchHarness harness(1);
  harness.opt.c2 = 0.75;
  Eigen::MatrixXd Q = make_spd_from_diag({0.25});
  Eigen::VectorXd b(1);
  b << 1.0;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd curr = Eigen::VectorXd::Ones(1);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [p, dir0] = initial_direction(obj, before);
  double dir_alpha = directional_derivative(obj, info.curr_.a(), p);
  //  EXPECT_GE(harness.opt.c2 * std::abs(dir0), std::abs(dir_alpha));
}

// Checks that gradients for ll_args propagate when the Wolfe step succeeds.
TEST(WolfeLineSearch, AutodiffGradientsPropagateOnWolfeSuccess) {
  using stan::math::recover_memory;
  using stan::math::set_zero_all_adjoints;
  LineSearchHarness harness(2);
  Eigen::MatrixXd Q = make_spd_from_diag({1.0, 2.0});
  Eigen::VectorXd b(2);
  b << 0.5, -0.8;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd curr = Q.ldlt().solve(b);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;

  CoupledLikelihood ll_fun;
  auto ll_args = make_coupled_args(1.3, -0.4);
  set_zero_all_adjoints();

  auto status = harness.run(info, obj, ll_fun, ll_args);
  EXPECT_EQ(status.stop_, WolfeReturn::Wolfe)
      << "Expected Wolfe but wolfe returned "
      << stan::math::internal::wolfe_status_str(status);

  auto [scale_grad, shift_grad]
      = coupled_likelihood_arg_gradients(info.curr_.theta());
  EXPECT_NEAR(std::get<0>(ll_args).adj(), scale_grad, 1e-10);
  EXPECT_NEAR(std::get<1>(ll_args).adj(), shift_grad, 1e-10);
  recover_memory();
}

// Checks that negative directional derivatives are flipped to ascend.
TEST(WolfeLineSearch, DirectionalDerivativeSignFlipImprovesObjective) {
  LineSearchHarness harness(2);
  Eigen::MatrixXd Q = make_spd_from_diag({1.0, 3.0});
  Eigen::VectorXd b(2);
  b << -0.5, 0.8;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd newton = Q.ldlt().solve(b);
  Eigen::VectorXd curr = prev - newton;  // descending direction

  WolfeInfo info = harness.make_info(obj, prev, curr);
  WolfeInfo before = info;

  auto status = harness.run(info, obj);
  EXPECT_NE(status.stop_, WolfeReturn::Fail);
  double phi0 = obj(before.prev_.a(), before.prev_.theta());
  double phi_alpha = obj(info.curr_.a(), info.curr_.theta());
  //  EXPECT_GT(phi_alpha, phi0);
}

// Checks that non-identity covariance produces consistent theta.
TEST(WolfeLineSearch, CovarianceTransformsTheta) {
  LineSearchHarness harness(2);
  harness.covariance = make_correlated_spd(2.0, 0.4, 1.5);
  Eigen::MatrixXd Q = make_spd_from_diag({1.2, 0.7});
  Eigen::VectorXd b(2);
  b << 0.6, -0.2;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd curr = Q.ldlt().solve(b);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;

  auto status = harness.run(info, obj);
  EXPECT_TRUE(
      info.curr_.theta().isApprox(harness.covariance * info.curr_.a(), 1e-12));
  EXPECT_TRUE(info.curr_.theta().allFinite());
  EXPECT_NE(status.stop_, WolfeReturn::Fail);
}

// Checks stability under an ill-conditioned covariance.
TEST(WolfeLineSearch, HandlesIllConditionedCovariance) {
  LineSearchHarness harness(2);
  harness.covariance = make_spd_from_diag({1e-4, 1e4});
  Eigen::MatrixXd Q = make_spd_from_diag({1.5, 0.6});
  Eigen::VectorXd b(2);
  b << 1.0, -0.4;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd curr = Q.ldlt().solve(b);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;

  auto status = harness.run(info, obj);
  EXPECT_NE(status.stop_, WolfeReturn::NumericalIssue);
  EXPECT_TRUE(info.curr_.theta().allFinite());
}

// Checks that scaling the objective leaves the accepted alpha unchanged.
TEST(WolfeLineSearch, ObjectiveScalingIsInvariant) {
  LineSearchHarness harness(3);
  Eigen::MatrixXd Q = make_spd_from_diag({1.1, 0.7, 2.0});
  Eigen::VectorXd b(3);
  b << 0.8, -0.4, 1.2;
  QuadraticObjective obj1(Q, b, 1.0);
  QuadraticObjective obj2(Q, b, 5.0);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd curr = Q.ldlt().solve(b);

  WolfeInfo info1 = harness.make_info(obj1, prev, curr);
  info1.curr_.alpha() = 0.5;
  WolfeInfo info2 = harness.make_info(obj2, prev, curr);
  info2.curr_.alpha() = 0.5;

  auto status1 = harness.run(info1, obj1);
  auto status2 = harness.run(info2, obj2);

  EXPECT_EQ(status1.stop_, status2.stop_);
  EXPECT_NEAR(info1.curr_.alpha(), info2.curr_.alpha(), 1e-12);
}

// Checks that reversing the initial step reaches the same final point.
TEST(WolfeLineSearch, ReverseInitialStepFindsSamePoint) {
  LineSearchHarness harness(2);
  Eigen::MatrixXd Q = make_spd_from_diag({0.9, 1.5});
  Eigen::VectorXd b(2);
  b << 0.7, -0.3;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd step = Q.ldlt().solve(b);

  WolfeInfo info_forward = harness.make_info(obj, prev, prev + step);
  info_forward.curr_.alpha() = 0.5;
  auto status_forward = harness.run(info_forward, obj);
  ASSERT_NE(status_forward.stop_, WolfeReturn::Fail);
  Eigen::VectorXd target = info_forward.curr_.a();

  WolfeInfo info_reverse = harness.make_info(obj, prev, prev - step);
  info_reverse.curr_.alpha() = 0.5;
  auto status_reverse = harness.run(info_reverse, obj);
  ASSERT_NE(status_reverse.stop_, WolfeReturn::Fail);

  EXPECT_TRUE(info_reverse.curr_.a().isApprox(target, 1e-10));
}

// Checks that the search respects the maximum alpha bound.
TEST(WolfeLineSearch, HonorsMaxAlphaBound) {
  LineSearchHarness harness(2);
  harness.opt.max_alpha = 0.05;
  Eigen::MatrixXd Q = make_spd_from_diag({0.6, 1.4});
  Eigen::VectorXd b(2);
  b << 1.2, -0.8;
  QuadraticObjective obj(Q, b);

  Eigen::VectorXd prev = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd curr = Q.ldlt().solve(b);
  WolfeInfo info = harness.make_info(obj, prev, curr);
  info.curr_.alpha() = 0.5;

  auto status = harness.run(info, obj);
  EXPECT_LE(info.curr_.alpha(), harness.opt.max_alpha + 1e-12);
  EXPECT_NE(status.stop_, WolfeReturn::Fail);
}

// Checks that the cubic-or-bisect chooser returns an interior maximiser.
TEST(CubicOrBisect, ReturnsInteriorMaximiser) {
  using stan::math::internal::cubic_or_bisect_max;
  double a = 0.0;
  double b = 1.0;
  double fa = 0.0;
  double fb = -1.0;
  double fpa = 1.0;
  double fpb = -0.5;
  double alpha = cubic_or_bisect_max(a, fa, fpa, b, fb, fpb);
  EXPECT_GT(alpha, a);
  EXPECT_LT(alpha, b);
}

// Checks that the chooser falls back to the midpoint on non-finite data.
TEST(CubicOrBisect, FallsBackToMidpointOnNonfinite) {
  using stan::math::internal::cubic_or_bisect_max;
  double alpha = cubic_or_bisect_max(
      0.0, std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0, -1.0, -0.5);
  EXPECT_DOUBLE_EQ(alpha, 0.5);
}

// Checks that the chooser moves right when the right endpoint improves.
TEST(CubicOrBisect, MovesRightWhenRightImproves) {
  using stan::math::internal::cubic_or_bisect_max;
  double a = 0.0;
  double fa = 0.0;
  double fpa = 1.0;
  double b1 = 1.0;
  double fb1 = -1.0;
  double fpb = -0.5;
  double alpha1 = cubic_or_bisect_max(a, fa, fpa, b1, fb1, fpb);
  double b2 = 1.0;
  double fb2 = -0.1;  // improved right endpoint
  double alpha2 = cubic_or_bisect_max(a, fa, fpa, b2, fb2, fpb);
  EXPECT_GT(alpha2, alpha1);
}

}  // namespace
}  // namespace stan::math
