#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// Tests for integrate_1d_gauss_kronrod. Mirrors integrate_1d_test.cpp, with
// the following differences:
//   * functors that depended on the xc argument have been rewritten to use
//     the explicit distance-to-boundary expression instead, because
//     Gauss-Kronrod does not produce xc (it is always NaN here);
//   * the f11 xc==NaN sentinel test is omitted (xc is unconditionally NaN
//     under Gauss-Kronrod, so the original semantics do not apply);
//
// Note on the divide of labour vs integrate_1d:
//   - integrate_1d (tanh_sinh / exp_sinh / sinh_sinh, double-exponential
//     quadrature) excels at integrals with algebraic or logarithmic
//     endpoint singularities (e.g. 1/sqrt(x) near x=0, 1/sqrt(1-x) near
//     x=1, beta-type integrands x^{a-1}(1-x)^{b-1} with small a,b).
//   - Gauss-Kronrod has no endpoint transform and fails on those cases
//     unless the user pre-splits the interval; in exchange, it is faster
//     and more accurate on smooth integrands and handles modest
//     oscillation via adaptive bisection. The test cases below are
//     restricted to integrands where Gauss-Kronrod is competitive; the
//     endpoint-singular cases from the integrate_1d test suite are
//     deliberately omitted here.

namespace integrate_1d_gk_test {

std::ostringstream *msgs = nullptr;

struct f1 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(-x) / sqrt(x);
  }
};

// Original f2 used xc near x=1; rewritten with explicit (1 - x).
struct f2 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    if (x <= 0.5) {
      return sqrt(x) / sqrt(1 - x * x);
    } else {
      return sqrt(x / ((x + 1) * (1 - x)));
    }
  }
};

struct f3 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(-x);
  }
};

struct f4 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(x) + theta[0];
  }
};

struct f5 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(x) + pow(theta[0], 2) + pow(theta[1], 3);
  }
};

struct f6 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(x) + pow(x_i[0], 2) + pow(theta[0], 4) + 3 * theta[1];
  }
};

struct f7 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(x) + pow(x_r[0], 2) + pow(x_r[1], 5) + 3 * x_r[2];
  }
};

struct f8 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(-pow(x - theta[0], x_i[0]) / pow(x_r[0], x_i[0]));
  }
};

struct f9 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return 1.0 / (1.0 + pow(x, x_i[0]) / theta[0]);
  }
};

// Original f10 used xc on the right half; rewritten with explicit (1 - x).
struct f10 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return pow(x, theta[0] - 1.0) * pow(1 - x, theta[1] - 1.0);
  }
};

struct f12 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    T1 out = stan::math::modified_bessel_second_kind(0, x);
    if (out > 0)
      return 2 * x * out;
    return out;
  }
};

struct f13 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    T1 out = stan::math::modified_bessel_second_kind(0, x);
    if (out > 0)
      return 2 * x * stan::math::square(out);
    return out;
  }
};

// Original f14 used xc near x=1; rewritten with explicit (1 - x).
struct f14 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return exp(x) * stan::math::inv_sqrt(1 - x);
  }
};

struct f16 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<double> &x_r,
      const std::vector<int> &x_i) const {
    return x * sin(x) / (1 + stan::math::square(cos(x)));
  }
};

struct f17 {
  inline double operator()(const double &x, const double &xc,
                           std::ostream *msgs, const std::vector<double> &theta,
                           const std::vector<double> &x_r,
                           const std::vector<int> &x_i) const {
    double mu = theta[0];
    double sigma = theta[1];
    return 1.0 / (sqrt(2.0 * stan::math::pi()) * sigma)
           * std::exp(-0.5 * ((x - mu) / sigma) * ((x - mu) / sigma));
  }
};

/*
 * test_integration is a helper that integrates the provided function and
 * checks the computed value against val. It also exercises the flipped
 * domain (-b, -a) by negating x in the user functor.
 */
template <typename F>
inline void test_integration(const F &f, double a, double b,
                             std::vector<double> thetas,
                             const std::vector<double> &x_r,
                             const std::vector<int> &x_i, double val) {
  using stan::math::integrate_1d_gauss_kronrod;
  using stan::math::integrate_1d_gauss_kronrod_tol;

  std::vector<double> tolerances = {1e-4, 1e-6, 1e-8};

  for (auto tolerance : tolerances) {
    EXPECT_LE(std::abs(integrate_1d_gauss_kronrod_tol(
                           f, a, b, tolerance, 0.0, 15, msgs, thetas, x_r, x_i)
                       - val),
              tolerance);
    // Flip the domain of integration and check that the integral matches
    auto flipped
        = [&](const double &x, const double &xc, std::ostream *msgs,
              const std::vector<double> &theta, const std::vector<double> &x_r,
              const std::vector<int> &x_i) {
            return f(-x, -xc, msgs, theta, x_r, x_i);
          };
    EXPECT_LE(
        std::abs(integrate_1d_gauss_kronrod_tol(flipped, -b, -a, tolerance, 0.0,
                                                15, msgs, thetas, x_r, x_i)
                 - val),
        tolerance);
  }
}

}  // namespace integrate_1d_gk_test

TEST(StanMath_integrate_1d_gk_prim, TestThrows) {
  // Left limit of integration must be less than or equal to right limit
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f4{}, 1.0, 0.0, 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.5},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
  // NaN limits not okay
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f4{}, 0.0,
                   std::numeric_limits<double>::quiet_NaN(), 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.5},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_gauss_kronrod_tol(
          integrate_1d_gk_test::f4{}, std::numeric_limits<double>::quiet_NaN(),
          0.0, 1e-6, 0.0, 15, integrate_1d_gk_test::msgs,
          std::vector<double>{0.5}, std::vector<double>{}, std::vector<int>{}),
      std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_gauss_kronrod_tol(
          integrate_1d_gk_test::f4{}, std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN(), 1e-6, 0.0, 15,
          integrate_1d_gk_test::msgs, std::vector<double>{0.5},
          std::vector<double>{}, std::vector<int>{}),
      std::domain_error);
  // Two of the same inf limits not okay
  EXPECT_THROW(
      stan::math::integrate_1d_gauss_kronrod_tol(
          integrate_1d_gk_test::f4{}, -std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity(), 1e-6, 0.0, 15,
          integrate_1d_gk_test::msgs, std::vector<double>{0.5},
          std::vector<double>{}, std::vector<int>{}),
      std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_gauss_kronrod_tol(
          integrate_1d_gk_test::f4{}, std::numeric_limits<double>::infinity(),
          std::numeric_limits<double>::infinity(), 1e-6, 0.0, 15,
          integrate_1d_gk_test::msgs, std::vector<double>{0.5},
          std::vector<double>{}, std::vector<int>{}),
      std::domain_error);
  // Negative max_depth not okay
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f4{}, 0.0, 1.0, 1e-6, 0.0, -1,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.5},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
  // Negative absolute_tolerance not okay
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f4{}, 0.0, 1.0, 1e-6, -1e-3, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.5},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
}

TEST(StanMath_integrate_1d_gk_prim, test_integer_arguments) {
  // Use a smooth integrand for the integer-bounds smoke test; f4 is exp(x)+c
  // and integrates cleanly under Gauss-Kronrod.
  EXPECT_NO_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0, 1, 1e-6, 0.0, 15,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{}));
  EXPECT_NO_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0.0, 1, 1e-6, 0.0, 15,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{}));
  EXPECT_NO_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0, 1.0, 1e-6, 0.0, 15,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{}));
}

TEST(StanMath_integrate_1d_gk_prim, test1_smooth) {
  // Zero-crossing integral + limit at infinity (smooth exponential decay)
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f3{}, -2.0, std::numeric_limits<double>::infinity(),
      {}, {}, {}, 7.38905609893065);
  // Easy integrals
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f4{}, 0.2, 0.7, {0.5}, std::vector<double>{},
      std::vector<int>{}, 1.0423499493102901);
  integrate_1d_gk_test::test_integration(integrate_1d_gk_test::f5{}, -0.2, 0.7,
                                         {0.4, 0.4}, std::vector<double>{},
                                         std::vector<int>{}, 1.396621954392482);
  // Zero-length intervals
  integrate_1d_gk_test::test_integration(integrate_1d_gk_test::f4{}, 0.0, 0.0,
                                         {0.5}, std::vector<double>{},
                                         std::vector<int>{}, 0.0);
  integrate_1d_gk_test::test_integration(integrate_1d_gk_test::f5{}, 1.0, 1.0,
                                         {0.4, 0.4}, std::vector<double>{},
                                         std::vector<int>{}, 0.0);
  // Test x_i
  integrate_1d_gk_test::test_integration(integrate_1d_gk_test::f6{}, -0.2, 2.9,
                                         {6.0, 5.1}, {}, {4},
                                         4131.985414616364);
  // Test x_r
  integrate_1d_gk_test::test_integration(integrate_1d_gk_test::f7{}, -0.2, 2.9,
                                         {}, {4.0, 6.0, 5.1}, {},
                                         24219.985414616367);
  // Both limits at infinity + test x_r/x_i (smooth Gaussian-shaped)
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f8{}, -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), {5.0}, {1.7}, {2},
      3.013171546539377);
  // Both limits at infinity + test x_i (smooth rational on R)
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f9{}, -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), {1.3}, {}, {4},
      2.372032924895055);
  // Smooth oscillation
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f16{}, 0.0, stan::math::pi(), {}, {}, {},
      stan::math::square(stan::math::pi()) / 4);
  // Bounds working right (Gaussian PDF tail integral)
  integrate_1d_gk_test::test_integration(
      integrate_1d_gk_test::f17{}, -std::numeric_limits<double>::infinity(),
      -1.5, {0.0, 1.0}, {}, {}, 0.066807201268858071);
}

// Demonstrate the known weakness of Gauss-Kronrod: integrands with algebraic
// or logarithmic endpoint singularities (1/sqrt(x), 1/sqrt(1-x), beta-type
// densities with small parameters, ...). Boost's gauss_kronrod has no
// endpoint transform; without user-driven interval splitting it either
// converges very slowly or signals failure via the error estimate. The
// existing integrate_1d (tanh_sinh/exp_sinh) handles these cases natively
// and remains the preferred choice for them. This test documents the
// behaviour so future maintainers do not mistake it for a regression.
TEST(StanMath_integrate_1d_gk_prim, endpoint_singularity_throws) {
  // 1/sqrt(x) at x = 0 (f1)
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f1{}, 0.0,
                   std::numeric_limits<double>::infinity(), 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>(),
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
  // 1/sqrt(1-x*x) at x = 1 (f2)
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f2{}, 0.0, 1.0, 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>(),
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
  // beta integrand with small shape parameters (f10, a=b=0.1)
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f10{}, 0.0, 1.0, 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.1, 0.1},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);
}

TEST(StanMath_integrate_1d_gk_prim, max_depth_argument) {
  // Smoke test: explicit max_depth is accepted and produces a sensible result.
  // Argument order is (rel_tol, abs_tol, max_depth).
  double Q = stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0.2, 0.7, 1e-8, 0.0, 20,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{});
  EXPECT_NEAR(Q, 1.0423499493102901, 1e-8);
}

// Demonstrate that an explicit absolute_tolerance can suppress the
// convergence throw on integrands where the strict relative-tolerance
// test fails. We reuse f10 (beta integrand x^{a-1}(1-x)^{b-1}) with
// small shape parameters: this has algebraic endpoint singularities
// that integrate_1d_gauss_kronrod cannot resolve to the requested
// relative tolerance (covered by endpoint_singularity_throws above).
// Setting abs_tol large enough that
//   max(rel_tol * L1, abs_tol) >= reported_error
// lets the user accept the (possibly imprecise) estimate without an
// exception, matching the QUADPACK convention of mixed
// relative/absolute convergence.
TEST(StanMath_integrate_1d_gk_prim, abs_tol_suppresses_throw) {
  // Sanity: with abs_tol = 0 (default) the call throws (this is the
  // same case as endpoint_singularity_throws.f10 above).
  EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                   integrate_1d_gk_test::f10{}, 0.0, 1.0, 1e-6, 0.0, 15,
                   integrate_1d_gk_test::msgs, std::vector<double>{0.1, 0.1},
                   std::vector<double>{}, std::vector<int>{}),
               std::domain_error);

  // With a very generous abs_tol the convergence threshold is
  // satisfied and the integral is returned. The endpoint singularity
  // x^{-0.9}*(1-x)^{-0.9} makes the quadrature evaluate the integrand
  // at values approaching 1e9 near x=0, so the reported error estimate
  // is also large in absolute terms (~5e4 here); abs_tol = 1e6 is
  // safely above it. The true value of B(0.1, 0.1) is ~19.7, so even
  // an imprecise estimate should be in the right ballpark.
  double Q = 0.0;
  EXPECT_NO_THROW(Q = stan::math::integrate_1d_gauss_kronrod_tol(
                      integrate_1d_gk_test::f10{}, 0.0, 1.0, 1e-6, 1e6, 15,
                      integrate_1d_gk_test::msgs, std::vector<double>{0.1, 0.1},
                      std::vector<double>{}, std::vector<int>{}));
  EXPECT_GT(Q, 1.0);
  EXPECT_LT(Q, 1000.0);
}

TEST(StanMath_integrate_1d_gk_prim, abs_tol_argument_smoke) {
  // Smoke test: explicit abs_tol on a well-converged integrand does
  // not change the result. Argument order is (rel_tol, abs_tol).
  double Q0 = stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0.2, 0.7, 1e-8, 0.0, 15,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{});
  double Q1 = stan::math::integrate_1d_gauss_kronrod_tol(
      integrate_1d_gk_test::f4{}, 0.2, 0.7, 1e-8, 1e-12, 15,
      integrate_1d_gk_test::msgs, std::vector<double>{0.5},
      std::vector<double>{}, std::vector<int>{});
  EXPECT_NEAR(Q0, 1.0423499493102901, 1e-8);
  EXPECT_NEAR(Q1, Q0, 1e-12);
}

// ---------------------------------------------------------------------------
// absolute_tolerance during refinement
//
// absolute_tolerance is applied in two places, in the same units: as a floor
// on refinement (a panel already below the floor is not bisected) and as the
// floor on the convergence test. The tests below pin the three properties
// that make that safe:
//
//   1. equivalence  - abs_tol == 0 is bit-for-bit Boost;
//   2. continuity   - a negligible abs_tol is a no-op, not a mode switch;
//   3. monotonicity - raising abs_tol can only reduce work.
//
// Property 3 is the one that bounds the cost. An implementation that instead
// refines against a *global* error budget satisfies 1 and 2 but violates 3
// catastrophically: it drains the whole bisection tree whenever the tolerance
// is unreachable, which is exactly the regime abs_tol exists to serve.
// ---------------------------------------------------------------------------

// abs_tol == 0 must reproduce Boost's gauss_kronrod::integrate exactly --
// not merely to within a tolerance, but bit-for-bit in the estimate, the
// error, the L1 norm and the number of integrand evaluations. This is what
// licenses integrate_gk to route every call through the local
// implementation rather than dispatching to Boost when abs_tol == 0.
TEST(StanMath_integrate_1d_gk_prim,
     matches_boost_bit_for_bit_when_abs_tol_zero) {
  const double inf = std::numeric_limits<double>::infinity();

  auto expect_bit_exact = [](const char *name, auto f, double a, double b,
                             unsigned int depth, double rel_tol) {
    double boost_error = 0.0, boost_l1 = 0.0;
    double stan_error = 0.0, stan_l1 = 0.0;
    int boost_evaluations = 0, stan_evaluations = 0;

    auto boost_integrand = [&boost_evaluations, &f](double x) {
      ++boost_evaluations;
      return f(x);
    };
    auto stan_integrand = [&stan_evaluations, &f](double x) {
      ++stan_evaluations;
      return f(x);
    };

    const double boost_result
        = boost::math::quadrature::gauss_kronrod<double, 21>::integrate(
            boost_integrand, a, b, depth, rel_tol, &boost_error, &boost_l1);
    const double stan_result = stan::math::internal::gauss_kronrod_21_integrate(
        stan_integrand, a, b, depth, rel_tol, 0.0, &stan_error, &stan_l1);

    EXPECT_EQ(boost_result, stan_result) << name << ": estimate";
    EXPECT_EQ(boost_error, stan_error) << name << ": error";
    EXPECT_EQ(boost_l1, stan_l1) << name << ": L1";
    EXPECT_EQ(boost_evaluations, stan_evaluations) << name << ": evaluations";
  };

  // Finite limits, spanning easy, oscillatory and singular integrands.
  expect_bit_exact(
      "smooth", [](double x) { return std::exp(x); }, 0.0, 1.0, 15, 1e-10);
  expect_bit_exact(
      "narrow interval", [](double x) { return std::sin(1 / (x + 0.01)); }, 0.2,
      0.7, 15, 1e-8);
  expect_bit_exact(
      "oscillatory", [](double x) { return std::sin(127 * x); }, 0.0, 1.0, 15,
      1e-12);
  expect_bit_exact(
      "negligible amplitude",
      [](double x) { return 1e-12 * std::sin(127 * x); }, 0.0, 1.0, 5, 1e-12);
  expect_bit_exact(
      "endpoint singularity", [](double x) { return std::pow(x, -0.9); },
      1e-300, 1.0, 15, 1e-12);
  expect_bit_exact(
      "peaked", [](double x) { return 1 / (1 + 1e4 * x * x); }, -1.0, 1.0, 15,
      1e-10);
  expect_bit_exact(
      "identically zero", [](double x) { return 0.0; }, 0.0, 1.0, 15, 1e-12);
  expect_bit_exact(
      "max_depth zero", [](double x) { return std::exp(x); }, 0.0, 1.0, 0,
      1e-14);

  // All three infinite-limit changes of variable.
  expect_bit_exact(
      "right infinite", [](double x) { return std::exp(-x); }, 0.0, inf, 15,
      1e-10);
  expect_bit_exact(
      "right infinite heavy tail", [](double x) { return 1 / (1 + x * x); },
      0.0, inf, 15, 1e-12);
  expect_bit_exact(
      "left infinite", [](double x) { return std::exp(x); }, -inf, 0.0, 15,
      1e-10);
  expect_bit_exact(
      "doubly infinite", [](double x) { return std::exp(-x * x); }, -inf, inf,
      15, 1e-10);
  expect_bit_exact(
      "doubly infinite heavy tail", [](double x) { return 1 / (1 + x * x); },
      -inf, inf, 15, 1e-12);
}

// An absolute tolerance far below the relative target must be a no-op: it
// cannot change the answer and cannot change the amount of work. Without
// this, abs_tol == 0 is a sentinel rather than a limit, and "pass a tiny
// abs_tol to be safe" silently selects different behaviour.
TEST(StanMath_integrate_1d_gk_prim, abs_tol_is_continuous_at_zero) {
  auto run = [](double absolute_tolerance, int *evaluations) {
    auto integrand = [evaluations](double x, double xc, std::ostream *msgs) {
      ++*evaluations;
      return std::exp(-x * x) * std::cos(30 * x);
    };
    return stan::math::integrate_1d_gauss_kronrod_tol(
        integrand, 0.0, 3.0, 1e-10, absolute_tolerance, 15,
        integrate_1d_gk_test::msgs);
  };

  int zero_evaluations = 0, tiny_evaluations = 0;
  const double zero_result = run(0.0, &zero_evaluations);
  const double tiny_result = run(1e-300, &tiny_evaluations);

  EXPECT_EQ(zero_result, tiny_result);
  EXPECT_EQ(zero_evaluations, tiny_evaluations);
}

// The motivating case: an integrand whose magnitude is so small that the
// relative-tolerance test degenerates into comparing accumulated round-off
// against itself. A positive abs_tol stops the pointless refinement, and the
// answer is unchanged to well within the tolerance the caller asked for.
TEST(StanMath_integrate_1d_gk_prim,
     positive_abs_tol_reduces_work_on_negligible_integrand) {
  constexpr double scale = 1e-12;
  constexpr double frequency = 127.0;
  constexpr double absolute_tolerance = 1e-14;
  const double expected = scale * (1.0 - std::cos(frequency)) / frequency;

  auto run = [](double abs_tol, int *evaluations) {
    auto integrand = [evaluations](double x, double xc, std::ostream *msgs) {
      ++*evaluations;
      return scale * std::sin(frequency * x);
    };
    return stan::math::integrate_1d_gauss_kronrod_tol(
        integrand, 0.0, 1.0, 1e-12, abs_tol, 5, integrate_1d_gk_test::msgs);
  };

  int relative_evaluations = 0, absolute_evaluations = 0;
  const double relative_result = run(0.0, &relative_evaluations);
  const double absolute_result = run(absolute_tolerance, &absolute_evaluations);

  // The contract: the answer is within the absolute tolerance requested.
  // Asserting anything tighter would be asserting an accident of how much
  // more accurate K21 happens to be than the caller asked for.
  EXPECT_NEAR(absolute_result, expected, absolute_tolerance);
  EXPECT_NEAR(relative_result, expected, absolute_tolerance);
  EXPECT_NEAR(absolute_result, relative_result, absolute_tolerance);

  // Both paths also happen to be near machine precision here; check that
  // relatively rather than against a hard-coded absolute epsilon.
  EXPECT_LT(std::abs(absolute_result - expected) / std::abs(expected), 1e-10);

  // The point of the exercise: materially less work. Stated as a ratio so
  // this does not pin the exact recursion counts of the quadrature.
  EXPECT_LT(absolute_evaluations, relative_evaluations);
  EXPECT_LT(2 * absolute_evaluations, relative_evaluations);
}

// Regression guard for unbounded refinement.
//
// x^{-0.9} has an endpoint singularity that Gauss-Kronrod cannot resolve, so
// no attainable tolerance is ever met and every panel looks "not yet good
// enough". Refinement must still be driven panel-by-panel, so the cost stays
// at the level of the abs_tol == 0 call (~1.7e3 evaluations at the default
// max_depth). An implementation that refines against a global budget instead
// exhausts the entire depth-15 tree here: ~1.4e6 evaluations, three orders of
// magnitude more, for an identical answer.
TEST(StanMath_integrate_1d_gk_prim,
     positive_abs_tol_bounds_work_on_unresolvable_integrand) {
  auto run = [](double absolute_tolerance, int *evaluations) {
    auto integrand = [evaluations](double x, double xc, std::ostream *msgs) {
      ++*evaluations;
      return std::pow(x, -0.9);
    };
    // Not resolvable to this tolerance, so the convergence test fails and
    // the call throws; the evaluation count is what is under test.
    EXPECT_THROW(stan::math::integrate_1d_gauss_kronrod_tol(
                     integrand, 1e-300, 1.0, 1e-12, absolute_tolerance, 15,
                     integrate_1d_gk_test::msgs),
                 std::domain_error);
  };

  int relative_evaluations = 0, absolute_evaluations = 0;
  run(0.0, &relative_evaluations);
  run(1e-14, &absolute_evaluations);

  // Monotonicity: a floor on refinement can only remove work, never add it.
  EXPECT_LE(absolute_evaluations, relative_evaluations);

  // Absolute backstop, three orders of magnitude below the global-budget
  // failure mode and one order above the actual cost.
  EXPECT_LT(absolute_evaluations, 50000);
}

// Raising abs_tol must never increase the work done, on any integrand.
TEST(StanMath_integrate_1d_gk_prim, work_is_monotone_in_abs_tol) {
  auto evaluations_for = [](double absolute_tolerance) {
    int evaluations = 0;
    auto integrand = [&evaluations](double x, double xc, std::ostream *msgs) {
      ++evaluations;
      return std::exp(-x * x) * std::cos(30 * x);
    };
    try {
      stan::math::integrate_1d_gauss_kronrod_tol(integrand, 0.0, 3.0, 1e-12,
                                                 absolute_tolerance, 12,
                                                 integrate_1d_gk_test::msgs);
    } catch (const std::domain_error &) {
      // Convergence failure is irrelevant here; only the cost is.
    }
    return evaluations;
  };

  int previous = evaluations_for(0.0);
  for (double absolute_tolerance : {1e-300, 1e-16, 1e-12, 1e-8, 1e-4, 1e-1}) {
    const int current = evaluations_for(absolute_tolerance);
    EXPECT_LE(current, previous) << "abs_tol = " << absolute_tolerance;
    previous = current;
  }
}

// A positive abs_tol must not paper over a genuinely unconverged result: if
// the error estimate still exceeds max(rel_tol * L1, abs_tol), it throws.
TEST(StanMath_integrate_1d_gk_prim, positive_abs_tol_still_throws_when_needed) {
  auto integrand = [](double x, double xc, std::ostream *msgs) {
    return std::pow(x, -0.9);
  };
  EXPECT_THROW(
      stan::math::integrate_1d_gauss_kronrod_tol(
          integrand, 1e-300, 1.0, 1e-12, 1e-14, 15, integrate_1d_gk_test::msgs),
      std::domain_error);
}

// max_depth = 0 means "one panel, no bisection"; a positive abs_tol must not
// disturb that (the refinement floor is never consulted).
TEST(StanMath_integrate_1d_gk_prim, positive_abs_tol_with_max_depth_zero) {
  auto integrand
      = [](double x, double xc, std::ostream *msgs) { return std::exp(x); };
  const double Q = stan::math::integrate_1d_gauss_kronrod_tol(
      integrand, 0.0, 1.0, 1e-10, 1e-12, 0, integrate_1d_gk_test::msgs);
  EXPECT_NEAR(Q, std::exp(1.0) - 1.0, 1e-12);
}

// A positive abs_tol must agree with the abs_tol == 0 answer under every
// change of variable, not just the finite one.
TEST(StanMath_integrate_1d_gk_prim, positive_abs_tol_domain_transformations) {
  constexpr double relative_tolerance = 1e-12;
  constexpr double absolute_tolerance = 1e-10;
  constexpr int max_depth = 8;
  const double infinity = std::numeric_limits<double>::infinity();

  auto check_integral = [&](const auto &integrand, double lower, double upper,
                            double expected) {
    const double legacy_result = stan::math::integrate_1d_gauss_kronrod_tol(
        integrand, lower, upper, relative_tolerance, 0.0, max_depth,
        integrate_1d_gk_test::msgs);
    const double absolute_result = stan::math::integrate_1d_gauss_kronrod_tol(
        integrand, lower, upper, relative_tolerance, absolute_tolerance,
        max_depth, integrate_1d_gk_test::msgs);
    EXPECT_NEAR(absolute_result, expected, absolute_tolerance);
    EXPECT_NEAR(absolute_result, legacy_result, absolute_tolerance);
  };

  auto increasing_exponential
      = [](double x, double xc, std::ostream *msgs) { return std::exp(x); };
  auto decreasing_exponential
      = [](double x, double xc, std::ostream *msgs) { return std::exp(-x); };
  auto gaussian_kernel = [](double x, double xc, std::ostream *msgs) {
    return std::exp(-x * x);
  };

  check_integral(increasing_exponential, 0.0, 1.0, std::exp(1.0) - 1.0);
  check_integral(decreasing_exponential, 0.0, infinity, 1.0);
  check_integral(increasing_exponential, -infinity, 0.0, 1.0);
  check_integral(gaussian_kernel, -infinity, infinity,
                 std::sqrt(stan::math::pi()));
}
