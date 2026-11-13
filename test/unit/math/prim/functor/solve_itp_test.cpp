#include <stan/math/prim/functor/solve_itp.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <ostream>

namespace {

struct cubic_functor {
  mutable int evaluations = 0;

  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    ++evaluations;
    return x * x * x - x - 2.0;
  }
};

struct lambert_functor {
  mutable int evaluations = 0;

  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    ++evaluations;
    return x * exp(x) - 1.0;
  }
};

struct affine_functor {
  template <typename T, typename T_scale, typename T_offset>
  auto operator()(const T& x, std::ostream*, const T_scale& scale,
                  const T_offset& offset) const {
    return scale * x - offset;
  }
};

struct endpoint_functor {
  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    return x;
  }
};

struct no_bracket_functor {
  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    return x * x + 1.0;
  }
};

struct nonfinite_functor {
  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    if (x > 0.25 && x < 0.75) {
      return T(std::numeric_limits<double>::quiet_NaN());
    }
    return x - 0.5;
  }
};

struct steep_functor {
  mutable int evaluations = 0;

  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    ++evaluations;
    using std::atan;
    return atan(1e12 * (x - 0.123456789));
  }
};

struct final_width_roundoff_functor {
  mutable int evaluations = 0;

  template <typename T, typename... Args>
  auto operator()(const T& x, std::ostream*, const Args&...) const {
    ++evaluations;
    constexpr double root = 0.42425901218220274;
    constexpr double linear = 0.0004533163653989277;
    const auto z = x - root;
    return z * z * z + linear * z;
  }
};

}  // namespace

TEST(MathPrimFunctorSolveItp, BisectionIterationCount) {
  using stan::math::internal::itp_bisection_iterations;
  EXPECT_EQ(0, itp_bisection_iterations(1.0, 0.5));
  EXPECT_EQ(1, itp_bisection_iterations(1.0, 0.25));
  EXPECT_EQ(3, itp_bisection_iterations(8.0, 0.5));
  EXPECT_EQ(10, itp_bisection_iterations(1.0, 1.0 / 2048.0));
}

TEST(MathPrimFunctorSolveItp, PublishedCubicExample) {
  cubic_functor f;
  const double root
      = stan::math::solve_itp_tol(f, 1.0, 2.0, 5e-4, 0.1, 2.0, 1, 100, nullptr);

  EXPECT_NEAR(1.5213797068045676, root, 5e-4);
  // Two endpoint evaluations plus five ITP evaluations.
  EXPECT_EQ(7, f.evaluations);
}

TEST(MathPrimFunctorSolveItp, LambertExample) {
  lambert_functor f;
  const double root = stan::math::solve_itp_tol(f, -1.0, 1.0, 1e-10, 0.1, 2.0,
                                                1, 100, nullptr);

  EXPECT_NEAR(0.5671432904097838, root, 1e-10);
  EXPECT_LE(f.evaluations - 2, 8);
}

TEST(MathPrimFunctorSolveItp, VariadicArguments) {
  const double root
      = stan::math::solve_itp(affine_functor(), -10.0, 10.0, nullptr, 3.0, 2.0);
  EXPECT_NEAR(2.0 / 3.0, root, 1e-10);
}

TEST(MathPrimFunctorSolveItp, EndpointRoot) {
  EXPECT_EQ(0.0, stan::math::solve_itp(endpoint_functor(), 0.0, 2.0, nullptr));
}

TEST(MathPrimFunctorSolveItp, DefaultControlsOnSubnormalWidth) {
  const double denorm = std::numeric_limits<double>::denorm_min();
  EXPECT_EQ(
      0.0, stan::math::solve_itp(endpoint_functor(), -denorm, denorm, nullptr));
}

TEST(MathPrimFunctorSolveItp, MinmaxEvaluationBound) {
  steep_functor f;
  constexpr double lower = -1.0;
  constexpr double upper = 1.0;
  constexpr double x_tolerance = 1e-12;
  constexpr int64_t n0 = 0;
  const int64_t n_half = stan::math::internal::itp_bisection_iterations(
      upper - lower, x_tolerance);

  const double root = stan::math::solve_itp_tol(
      f, lower, upper, x_tolerance, 0.1, 2.0, n0, n_half + 1, nullptr);

  EXPECT_NEAR(0.123456789, root, x_tolerance);
  EXPECT_LE(f.evaluations, n_half + 3);
}

TEST(MathPrimFunctorSolveItp, RepairsOutwardRoundedFinalWidth) {
  final_width_roundoff_functor f;
  constexpr double x_tolerance = 1e-10;
  const int64_t exact_bound
      = stan::math::internal::itp_bisection_iterations(2.0, x_tolerance) + 1;

  const double root = stan::math::solve_itp_tol(
      f, -1.0, 1.0, x_tolerance, 0.1, 2.0, 1, exact_bound + 1, nullptr);

  EXPECT_NEAR(0.42425901218220274, root, x_tolerance);
  EXPECT_LE(f.evaluations, exact_bound + 3);
}

TEST(MathPrimFunctorSolveItp, RejectsNonBracket) {
  EXPECT_THROW(stan::math::solve_itp(no_bracket_functor(), -1.0, 1.0, nullptr),
               std::domain_error);
}

TEST(MathPrimFunctorSolveItp, RejectsInvalidControls) {
  cubic_functor f;
  EXPECT_THROW(
      stan::math::solve_itp_tol(f, 1.0, 2.0, 0.0, 0.1, 2.0, 1, 100, nullptr),
      std::domain_error);
  EXPECT_THROW(
      stan::math::solve_itp_tol(f, 1.0, 2.0, 1e-6, 0.0, 2.0, 1, 100, nullptr),
      std::domain_error);
  EXPECT_THROW(
      stan::math::solve_itp_tol(f, 1.0, 2.0, 1e-6, 0.1, 3.0, 1, 100, nullptr),
      std::domain_error);
  EXPECT_THROW(
      stan::math::solve_itp_tol(f, 1.0, 2.0, 1e-6, 0.1, 2.0, -1, 100, nullptr),
      std::domain_error);
}

TEST(MathPrimFunctorSolveItp, RejectsTooFewSteps) {
  cubic_functor f;
  EXPECT_THROW(
      stan::math::solve_itp_tol(f, 1.0, 2.0, 1e-12, 0.1, 2.0, 1, 2, nullptr),
      std::domain_error);
}

TEST(MathPrimFunctorSolveItp, RejectsNonfiniteFunctionValue) {
  EXPECT_THROW(stan::math::solve_itp_tol(nonfinite_functor(), 0.0, 1.0, 1e-8,
                                         0.2, 2.0, 1, 100, nullptr),
               std::domain_error);
}
