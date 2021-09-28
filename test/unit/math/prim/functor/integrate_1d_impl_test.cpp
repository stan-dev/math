#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace integrate_1d_impl_test {

std::ostringstream *msgs = nullptr;

struct f1 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    return exp(-x) / sqrt(x);
  }
};

struct f2 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    if (x <= 0.5) {
      return sqrt(x) / sqrt(1 - x * x);
    } else {
      return sqrt(x / ((x + 1) * (xc)));
    }
  }
};

struct f3 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    return exp(-x);
  }
};

struct f4 {
  template <typename T1, typename T2>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs,
                       const std::vector<T2> &theta) const {
    return exp(x) + theta[0];
  }
};

struct f5 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta) const {
    return exp(x) + pow(theta[0], 2) + pow(theta[1], 3);
  }
};

struct f6 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta, const std::vector<int> &x_i) const {
    return exp(x) + pow(x_i[0], 2) + pow(theta[0], 4) + 3 * theta[1];
  }
};

struct f7 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs,
                       const std::vector<double> &x_r) const {
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
      const std::vector<T2> &theta, const std::vector<int> &x_i) const {
    return 1.0 / (1.0 + pow(x, x_i[0]) / theta[0]);
  }
};

struct f10 {
  template <typename T1, typename T2>
  inline stan::return_type_t<T1, T2> operator()(
      const T1 &x, const T1 &xc, std::ostream *msgs,
      const std::vector<T2> &theta) const {
    return pow(x, theta[0] - 1.0)
           * pow((x > 0.5) ? xc : (1 - x), theta[1] - 1.0);
  }
};

struct f11 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    return (std::isnan(xc)) ? xc : 0.0;
  }
};

struct f12 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    T1 out = stan::math::modified_bessel_second_kind(0, x);
    if (out > 0)
      return 2 * x * out;
    return out;
  }
};

struct f13 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    T1 out = stan::math::modified_bessel_second_kind(0, x);
    if (out > 0)
      return 2 * x * stan::math::square(out);
    return out;
  }
};

struct f14 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    return exp(x) * stan::math::inv_sqrt(x > 0.5 ? xc : 1 - x);
  }
};

struct f15 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    T1 x2 = x * x;
    T1 numer = x2 * log(x);
    T1 denom = x < 0.5 ? (x + 1) * (x - 1) : (x + 1) * -xc;
    denom *= x2 * x2 + 1;
    return numer / denom;
  }
};

struct f16 {
  template <typename T1>
  inline T1 operator()(const T1 &x, const T1 &xc, std::ostream *msgs) const {
    return x * sin(x) / (1 + stan::math::square(cos(x)));
  }
};

double lbaX_pdf(double X, double t, double A, double v, double s,
                std::ostream *pstream__) {
  double b_A_tv_ts;
  double b_tv_ts;
  double term_1;
  double term_2;
  double pdf;

  b_A_tv_ts = (((X - A) - (t * v)) / (t * s));
  b_tv_ts = ((X - (t * v)) / (t * s));
  term_1 = stan::math::Phi(b_A_tv_ts);
  term_2 = stan::math::Phi(b_tv_ts);
  pdf = ((1 / A) * (-term_1 + term_2));
  return pdf;
}

double lbaX_cdf(double X, double t, double A, double v, double s,
                std::ostream *pstream__) {
  double b_A_tv;
  double b_tv;
  double ts;
  double term_1;
  double term_2;
  double term_3;
  double term_4;
  double cdf;

  b_A_tv = ((X - A) - (t * v));
  b_tv = (X - (t * v));
  ts = (t * s);
  term_1 = (b_A_tv * stan::math::Phi((b_A_tv / ts)));
  term_2 = (b_tv * stan::math::Phi((b_tv / ts)));
  term_3 = (ts
            * stan::math::exp(
                  stan::math::normal_lpdf<false>((b_A_tv / ts), 0, 1)));
  term_4
      = (ts
         * stan::math::exp(stan::math::normal_lpdf<false>((b_tv / ts), 0, 1)));
  cdf = ((1 / A) * (((-term_1 + term_2) - term_3) + term_4));
  return cdf;
}

double rank_density(double x, double xc, const std::vector<double> &theta,
                    const std::vector<double> &x_r, const std::vector<int> &x_i,
                    std::ostream *pstream__) {
  double t = theta[0];
  double A = theta[1];
  double v1 = theta[2];
  double v2 = theta[3];
  double s = theta[4];
  double v = (lbaX_pdf(x, t, A, v1, s, pstream__)
              * lbaX_cdf(x, t, A, v2, s, pstream__));
  return v;
}

struct rank_density_functor__ {
  double operator()(double x, double xc, std::ostream *pstream__,
                    const std::vector<double> &theta,
                    const std::vector<double> &x_r,
                    const std::vector<int> &x_i) const {
    return rank_density(x, xc, theta, x_r, x_i, pstream__);
  }
};

double order(double down, double up, const std::vector<double> &theta,
             const std::vector<double> &x_r, std::ostream *pstream__) {
  std::vector<int> x_i;

  double v;

  v = stan::math::integrate_1d_impl(rank_density_functor__(), down, up, 1e-8,
                                    pstream__, theta, x_r, x_i);
  return v;
}
}  // namespace integrate_1d_impl_test
/*
 * test_integration is a helper function to make it easy to test the
 * integrate_1d_impl function.
 *
 * It takes in a callable function object, integration limits, parameters, real
 * and integer data. It integrates the provided function and compares the
 * computed integral against the provided integral (val).
 *
 * The prototype for f is:
 *   struct f10 {
 *     inline double operator()(const double& x, const double& xc, const
 * std::vector<double>& theta, const std::vector<double>& x_r, const
 * std::vector<int>& x_i, std::ostream& msgs) const {
 *     }
 *   };
 *
 * @tparam F Type of f
 * @param f a functor with signature given above
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param param parameters to be passed to f (should be
 * std::vector<stan::math::var>)
 * @param x_r real data to be passed to f (should be std::vector<double>)
 * @param x_i integer data to be passed to f (should be std::vector<int>)
 * @param val correct value of integral
 */
template <typename F, typename... Args>
void test_integration(const F &f, double a, double b, double val,
                      const Args &... args) {
  using stan::math::integrate_1d_impl;

  std::vector<double> tolerances = {1e-4, 1e-6, 1e-8};

  for (auto tolerance : tolerances) {
    EXPECT_LE(std::abs(integrate_1d_impl(f, a, b, tolerance,
                                         integrate_1d_impl_test::msgs, args...)
                       - val),
              tolerance);
    // Flip the domain of integration and check that the integral is working
    auto flipped
        = [&](double x, double xc, std::ostream *msgs, const auto &... args) {
            return f(-x, -xc, msgs, args...);
          };

    EXPECT_LE(std::abs(integrate_1d_impl(flipped, -b, -a, tolerance,
                                         integrate_1d_impl_test::msgs, args...)
                       - val),
              tolerance);
  }
}

TEST(StanMath_integrate_1d_impl_prim, TestThrows) {
  // Left limit of integration must be less than or equal to right limit
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{}, 1.0, 0.0, 0.0,
                                    integrate_1d_impl_test::msgs),
      std::domain_error);
  // NaN limits not okay
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{}, 0.0,
                                    std::numeric_limits<double>::quiet_NaN(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{},
                                    std::numeric_limits<double>::quiet_NaN(),
                                    0.0, 0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{},
                                    std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::quiet_NaN(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  // Two of the same inf limits not okay
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{},
                                    -std::numeric_limits<double>::infinity(),
                                    -std::numeric_limits<double>::infinity(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::domain_error);

  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{},
                                    std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::infinity(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  // xc should be nan if there are infinite limits
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f11{}, 0.0,
                                    std::numeric_limits<double>::infinity(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::runtime_error);
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f11{},
                                    std::numeric_limits<double>::infinity(),
                                    0.0, 0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  EXPECT_THROW(
      stan::math::integrate_1d_impl(integrate_1d_impl_test::f11{},
                                    std::numeric_limits<double>::infinity(),
                                    std::numeric_limits<double>::infinity(),
                                    0.0, integrate_1d_impl_test::msgs),
      std::domain_error);
  // But not otherwise
  EXPECT_NO_THROW(stan::math::integrate_1d_impl(integrate_1d_impl_test::f11{},
                                                0.0, 1.0, 0.0,
                                                integrate_1d_impl_test::msgs));
}

TEST(StanMath_integrate_1d_impl_prim, test_integer_arguments) {
  double v;
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{}, 0, 1, 0.0,
                                        integrate_1d_impl_test::msgs));
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{}, 0.0, 1,
                                        0.0, integrate_1d_impl_test::msgs));
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d_impl(integrate_1d_impl_test::f2{}, 0, 1.0,
                                        0.0, integrate_1d_impl_test::msgs));
}

TEST(StanMath_integrate_1d_impl_prim, test1) {
  // Tricky integral from Boost docs + limit at infinity
  test_integration(integrate_1d_impl_test::f1{}, 0.0,
                   std::numeric_limits<double>::infinity(), 1.772453850905516);
  // Tricky integral from Boost 1d integration docs
  test_integration(integrate_1d_impl_test::f2{}, 0.0, 1.0, 1.198140234735592);
  // Tricky integral from Boost 1d integration docs
  test_integration(integrate_1d_impl_test::f2{}, 0, 1, 1.198140234735592);
  // Zero crossing integral + limit at infinity
  test_integration(integrate_1d_impl_test::f3{}, -2.0,
                   std::numeric_limits<double>::infinity(), 7.38905609893065);
  // Easy integrals
  test_integration(integrate_1d_impl_test::f4{}, 0.2, 0.7, 1.0423499493102901,
                   std::vector<double>({0.5}));
  test_integration(integrate_1d_impl_test::f5{}, -0.2, 0.7, 1.396621954392482,
                   std::vector<double>({0.4, 0.4}));
  test_integration(integrate_1d_impl_test::f4{}, 0.0, 0.0, 0.0,
                   std::vector<double>({0.5}));
  test_integration(integrate_1d_impl_test::f5{}, 1.0, 1.0, 0.0,
                   std::vector<double>({0.4, 0.4}));
  // Test x_i
  test_integration(integrate_1d_impl_test::f6{}, -0.2, 2.9, 4131.985414616364,
                   std::vector<double>({6.0, 5.1}), std::vector<int>({4}));
  // Test x_r
  test_integration(integrate_1d_impl_test::f7{}, -0.2, 2.9, 24219.985414616367,
                   std::vector<double>({4.0, 6.0, 5.1}));
  // Both limits at infinity + test x_r/x_i
  test_integration(integrate_1d_impl_test::f8{},
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(), 3.013171546539377,
                   std::vector<double>({5.0}), std::vector<double>({1.7}),
                   std::vector<int>({2}));
  // Both limits at infinity + test x_i
  test_integration(integrate_1d_impl_test::f9{},
                   -std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(), 2.372032924895055,
                   std::vector<double>({1.3}), std::vector<int>({4}));
  // Various integrals of beta function
  test_integration(integrate_1d_impl_test::f10{}, 0.0, 1.0, 19.71463948905016,
                   std::vector<double>({0.1, 0.1}));
  test_integration(integrate_1d_impl_test::f10{}, 0.0, 1.0, 11.32308697521577,
                   std::vector<double>({0.1, 0.5}));
  test_integration(integrate_1d_impl_test::f10{}, 0.0, 1.0, 11.32308697521577,
                   std::vector<double>({0.5, 0.1}));
  test_integration(integrate_1d_impl_test::f10{}, 0.0, 1.0, 0.00952380952380952,
                   std::vector<double>({5.0, 3.0}));

  // Integrals from
  // http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf
  test_integration(integrate_1d_impl_test::f12{}, 0.0,
                   std::numeric_limits<double>::infinity(), 2.0);
  test_integration(integrate_1d_impl_test::f13{}, 0.0,
                   std::numeric_limits<double>::infinity(), 1.0);
  test_integration(integrate_1d_impl_test::f14{}, 0.0, 1.0,
                   exp(1) * sqrt(stan::math::pi()) * stan::math::erf(1.0));

  // Integrals from http://crd-legacy.lbl.gov/~dhbailey/dhbpapers/quadrature.pdf
  // works normally but not to tolerance when limits of integration are flipped
  //  test_integration(f15{}, 0.0, 1.0, {}, {}, {},
  //                   stan::math::square(stan::math::pi()) * (2 - sqrt(2.0)) /
  //                   32);
  test_integration(integrate_1d_impl_test::f16{}, 0.0, stan::math::pi(),
                   stan::math::square(stan::math::pi()) / 4);
}

TEST(StanMath_integrate_1d_impl_prim, TestTolerance) {
  std::ostringstream *msgs = nullptr;

  double t = 0.5;
  double b = 1.0;
  double A = 0.5;
  double v1 = 1.0;
  double v2 = 1.0;
  double s = 1.0;

  std::vector<double> theta = {t, A, v1, v2, s};
  std::vector<double> x_r;

  EXPECT_NO_THROW(integrate_1d_impl_test::order(-10, 0.67, theta, x_r, msgs));
}
