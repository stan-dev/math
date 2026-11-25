#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <sstream>
#include <string>
#include <tuple>

namespace {

struct IdentityCovariance2 {
  Eigen::MatrixXd operator()(std::ostream* /*msgs*/) const {
    return Eigen::MatrixXd::Identity(2, 2);
  }
};

// Concave in the first coordinate, mildly convex in the second: the negative
// Hessian of the log likelihood has a negative diagonal entry, so solver 1
// (diagonal Hessian-root Cholesky) rejects it, while the posterior curvature
// 1 - 0.2 stays positive and solver 2 converges.
struct PartlyConvexLikelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, std::ostream* /*msgs*/) const {
    return -0.5 * stan::math::square(theta(0) - 1.0)
           + 0.1 * stan::math::square(theta(1));
  }
};

double run(bool allow_fallthrough, std::ostream* msgs) {
  Eigen::VectorXd theta0 = Eigen::VectorXd::Zero(2);
  return stan::math::laplace_marginal_tol<false>(
      PartlyConvexLikelihood{}, std::tuple<>{}, 1, IdentityCovariance2{},
      std::tuple<>{},
      std::make_tuple(theta0, 1e-10, 100, 1, 100, allow_fallthrough ? 1 : 0),
      msgs);
}

// Strictly concave everywhere, so each of the three solvers converges on
// its own and none of them has to fall back to the next one.
struct ConcaveLikelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, std::ostream* /*msgs*/) const {
    return -0.5 * stan::math::square(theta(0) - 1.0)
           - 0.5 * stan::math::square(theta(1) + 1.0);
  }
};

// The solver that produced the estimate for the given likelihood.
template <typename LL>
int used_solver(const LL& ll, int solver, bool allow_fallthrough) {
  Eigen::VectorXd theta0 = Eigen::VectorXd::Zero(2);
  stan::math::laplace_options_user_supplied options(
      theta0, 1e-10, 100, 1, solver, 100, allow_fallthrough);
  return stan::math::internal::laplace_marginal_density_est(
             ll, std::tuple<>{},
             Eigen::MatrixXd(Eigen::MatrixXd::Identity(2, 2)), options, nullptr)
      .solver_used;
}

}  // namespace

TEST(LaplaceSolverFallback, FallbackIsLoggedToStream) {
  std::ostringstream msgs;
  const double result = run(true, &msgs);
  EXPECT_TRUE(std::isfinite(result));
  EXPECT_NE(msgs.str().find("solver fallback"), std::string::npos);
  EXPECT_NE(msgs.str().find("not positive definite"), std::string::npos);
}

// A null message stream must not turn an allowed fallback into an error.
TEST(LaplaceSolverFallback, FallbackWithNullStreamDoesNotThrow) {
  double result = std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(result = run(true, nullptr));
  EXPECT_TRUE(std::isfinite(result));
}

// With fallthrough disabled the solver failure is reported as such, also
// after an earlier fallback in the same process has already been logged.
TEST(LaplaceSolverFallback, DisabledFallthroughThrowsWithReason) {
  std::ostringstream msgs;
  try {
    run(false, &msgs);
    FAIL() << "expected std::domain_error";
  } catch (const std::domain_error& e) {
    const std::string what = e.what();
    EXPECT_NE(what.find("allow_fallthrough"), std::string::npos) << what;
    EXPECT_NE(what.find("not positive definite"), std::string::npos) << what;
  }
  try {
    run(false, nullptr);
    FAIL() << "expected std::domain_error";
  } catch (const std::domain_error& e) {
    EXPECT_NE(std::string(e.what()).find("allow_fallthrough"),
              std::string::npos)
        << e.what();
  }
}

// The requested solver runs whether or not the fallthrough is allowed; the
// fallthrough only decides what happens after that solver has failed.
TEST(LaplaceSolverFallback, RequestedSolverRunsWithFallthroughAllowed) {
  for (int solver = 1; solver <= 3; ++solver) {
    EXPECT_EQ(solver, used_solver(ConcaveLikelihood{}, solver, true))
        << "solver " << solver;
    EXPECT_EQ(solver, used_solver(ConcaveLikelihood{}, solver, false))
        << "solver " << solver;
  }
  // Guard: after solver 1 fails, the fallthrough lands on solver 2.
  EXPECT_EQ(2, used_solver(PartlyConvexLikelihood{}, 1, true));
}
