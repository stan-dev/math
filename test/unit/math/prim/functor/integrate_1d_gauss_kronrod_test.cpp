#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
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

}  // namespace integrate_1d_gk_test

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
                           f, a, b, tolerance, 0.0, 15,
                           integrate_1d_gk_test::msgs, thetas, x_r, x_i)
                       - val),
              tolerance);
    // Flip the domain of integration and check that the integral matches
    auto flipped
        = [&](const double &x, const double &xc, std::ostream *msgs,
              const std::vector<double> &theta, const std::vector<double> &x_r,
              const std::vector<int> &x_i) {
            return f(-x, -xc, msgs, theta, x_r, x_i);
          };
    EXPECT_LE(std::abs(integrate_1d_gauss_kronrod_tol(
                           flipped, -b, -a, tolerance, 0.0, 15,
                           integrate_1d_gk_test::msgs, thetas, x_r, x_i)
                       - val),
              tolerance);
  }
}

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
  test_integration(integrate_1d_gk_test::f3{}, -2.0,
                   std::numeric_limits<double>::infinity(), {}, {}, {},
                   7.38905609893065);
  // Easy integrals
  test_integration(integrate_1d_gk_test::f4{}, 0.2, 0.7, {0.5},
                   std::vector<double>{}, std::vector<int>{},
                   1.0423499493102901);
  test_integration(integrate_1d_gk_test::f5{}, -0.2, 0.7, {0.4, 0.4},
                   std::vector<double>{}, std::vector<int>{},
                   1.396621954392482);
  // Zero-length intervals
  test_integration(integrate_1d_gk_test::f4{}, 0.0, 0.0, {0.5},
                   std::vector<double>{}, std::vector<int>{}, 0.0);
  test_integration(integrate_1d_gk_test::f5{}, 1.0, 1.0, {0.4, 0.4},
                   std::vector<double>{}, std::vector<int>{}, 0.0);
  // Test x_i
  test_integration(integrate_1d_gk_test::f6{}, -0.2, 2.9, {6.0, 5.1}, {}, {4},
                   4131.985414616364);
  // Test x_r
  test_integration(integrate_1d_gk_test::f7{}, -0.2, 2.9, {}, {4.0, 6.0, 5.1},
                   {}, 24219.985414616367);
  // Both limits at infinity + test x_r/x_i (smooth Gaussian-shaped)
  test_integration(integrate_1d_gk_test::f8{},
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(), {5.0}, {1.7}, {2},
                   3.013171546539377);
  // Both limits at infinity + test x_i (smooth rational on R)
  test_integration(integrate_1d_gk_test::f9{},
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(), {1.3}, {}, {4},
                   2.372032924895055);
  // Smooth oscillation
  test_integration(integrate_1d_gk_test::f16{}, 0.0, stan::math::pi(), {}, {},
                   {}, stan::math::square(stan::math::pi()) / 4);
  // Bounds working right (Gaussian PDF tail integral)
  test_integration(integrate_1d_gk_test::f17{},
                   -std::numeric_limits<double>::infinity(), -1.5, {0.0, 1.0},
                   {}, {}, 0.066807201268858071);
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
// lets the user accept Boost's (possibly imprecise) estimate without
// an exception, matching the QUADPACK convention of mixed
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
  // x^{-0.9}*(1-x)^{-0.9} makes Boost evaluate the integrand at
  // values approaching 1e9 near x=0, so the reported error estimate
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
