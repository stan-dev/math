#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>

#include <cmath>
#include <sstream>
#include <tuple>

namespace stan::math {
namespace {

struct IdentityCovariance {
  template <typename Stream>
  Eigen::MatrixXd operator()(Stream* /*msgs*/) const {
    return Eigen::MatrixXd::Identity(1, 1);
  }
};

struct QuarticLikelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, std::ostream* /*msgs*/) const {
    const auto& x = theta(0);
    const auto x_sq = stan::math::square(x);
    return 2.0 * x - 0.5 * x_sq - 0.5 * stan::math::square(x_sq);
  }
};

struct TinyQuarticLikelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, std::ostream* /*msgs*/) const {
    return 1e-8 * QuarticLikelihood{}(theta, nullptr);
  }
};

struct StubNewtonSolver {
  double proposal_a;

  template <typename NewtonStateT, typename LLFun, typename LLTupleArgs,
            typename CovarMat>
  void solve_step(NewtonStateT& state, const LLFun& /*ll_fun*/,
                  const LLTupleArgs& /*ll_args*/,
                  const CovarMat& /*covariance*/, int /*hessian_block_size*/,
                  std::ostream* /*msgs*/) const {
    state.proposal_step().a()(0) = proposal_a;
  }

  double compute_log_determinant() const { return 0.0; }

  template <typename NewtonStateT>
  double build_result(NewtonStateT& state, double /*log_det*/) const {
    return state.prev().a()(0);
  }
};

template <typename Likelihood>
double run_laplace(const Likelihood& ll_fun, double theta0_value,
                   double tolerance, int max_num_steps,
                   int max_steps_line_search, std::ostream* msgs) {
  Eigen::VectorXd theta0(1);
  theta0 << theta0_value;
  return stan::math::laplace_marginal_tol<false>(
      ll_fun, std::tuple<>{}, 1, IdentityCovariance{}, std::tuple<>{},
      std::make_tuple(theta0, tolerance, max_num_steps, 1,
                      max_steps_line_search, true),
      msgs);
}

TEST(LaplaceMarginalDensityEstimator, PublicLineSearchMatchesDirectStep) {
  std::ostringstream no_search_msgs;
  std::ostringstream wolfe_msgs;

  const double no_search
      = run_laplace(QuarticLikelihood{}, 2.0, 1e-12, 50, 0, &no_search_msgs);
  const double with_wolfe
      = run_laplace(QuarticLikelihood{}, 2.0, 1e-12, 50, 1000, &wolfe_msgs);

  EXPECT_TRUE(std::isfinite(no_search));
  EXPECT_TRUE(std::isfinite(with_wolfe));
  EXPECT_NEAR(no_search, with_wolfe, 1e-8);
}

TEST(LaplaceMarginalDensityEstimator, AbsoluteObjectiveToleranceStopsNearZero) {
  std::ostringstream msgs;

  const double result
      = run_laplace(TinyQuarticLikelihood{}, 0.0, 1e-8, 6, 1000, &msgs);

  EXPECT_TRUE(std::isfinite(result));
  EXPECT_EQ(msgs.str().find("max number of iterations"), std::string::npos);
}

TEST(LaplaceMarginalDensityEstimator,
     InvalidCachedProposalDoesNotTriggerArmijoFallback) {
  Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(1, 1);
  Eigen::VectorXd theta0 = Eigen::VectorXd::Zero(1);
  auto obj_fun = [](const auto& /*a*/, const auto& /*theta*/) { return -1.0; };
  auto theta_grad_f
      = [](const auto& theta) { return Eigen::VectorXd::Zero(theta.size()); };
  internal::NewtonState state(1, obj_fun, theta_grad_f, covariance, theta0);
  laplace_options_base options;
  options.hessian_block_size = 1;
  options.max_num_steps = 1;
  options.tolerance = 1e-12;
  options.line_search.max_iterations = 5;
  options.line_search.min_alpha = 1e-8;

  StubNewtonSolver solver{5.0};
  Eigen::Index step_iter = 1;
  auto failing_update = [min_alpha = options.line_search.min_alpha](
                            auto& /*proposal*/, auto&& /*curr*/,
                            auto&& /*prev*/, auto& eval_in, auto&& /*p*/) {
    eval_in.alpha() = 0.5 * min_alpha;
    return false;
  };
  auto unused_ll
      = [](const auto& /*theta*/, std::ostream* /*msgs*/) { return 0.0; };

  const double result = internal::run_newton_loop(
      solver, state, options, step_iter, unused_ll, std::tuple<>{}, covariance,
      failing_update, nullptr);

  EXPECT_DOUBLE_EQ(result, 0.0);
  EXPECT_FALSE(state.wolfe_status.accept_);
  EXPECT_EQ(state.wolfe_status.stop_, internal::WolfeReturn::StepTooSmall);
}

}  // namespace
}  // namespace stan::math
