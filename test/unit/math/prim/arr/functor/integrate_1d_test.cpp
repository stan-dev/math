#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>

std::ostringstream msgs;

struct f1 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(-x) / sqrt(x);
  }
};

struct f2 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return sqrt(x / (1 - x * x));
  }
};
  
struct f3 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(-x);
  }
};
  
struct f4 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + theta[0];
  }
};

struct f5 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + pow(theta[0], 2) + pow(theta[1], 3);
  }
};

struct f6 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + pow(x_i[0], 2) + pow(theta[0], 4) + 3 * theta[1];
  }
};

struct f7 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + pow(x_r[0], 2) + pow(x_r[1], 5) + 3 * x_r[2];
  }
};

struct f8 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(-pow(x - theta[0], x_i[0]) / pow(x_r[0], x_i[0]));
  }
};

struct f9 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return 1.0 / (1.0 + pow(x, x_i[0]) / theta[0]);
  }
};

/*
 * test_integration is a helper function to make it easy to test the integrate_1d
 * function.
 *
 * It takes in a callable function object, integration limits, parameters, real and
 * integer data. It integrates the provided function and compares the computed integral
 * against the provided integral (val).
 *
 * The prototype for f is:
 *   struct f10 {
 *     inline double operator()(const double& x, const std::vector<double>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
 *     }
 *  };
 *
 * @tparam F Type of f
 * @param f a functor with signature given above
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param param parameters to be passed to f (should be std::vector<stan::math::var>)
 * @param x_r real data to be passed to f (should be std::vector<double>)
 * @param x_i integer data to be passed to f (should be std::vector<int>)
 * @param val correct value of integral
 * @param grad correct value of gradients
 */
template <typename F>
void test_integration(const F& f, double a, double b, std::vector<double> thetas, const std::vector<double>& x_r, const std::vector<int>& x_i, double val) {
  using stan::math::integrate_1d;
  EXPECT_FLOAT_EQ(integrate_1d(f, a, b, thetas, x_r, x_i, msgs), val);
  // Flip the domain of integration and check that the integral is working
  auto flipped = [&](const double& x, const std::vector<double>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) { return f(-x, theta, x_r, x_i, msgs); };
  EXPECT_FLOAT_EQ(integrate_1d(flipped, -b, -a, thetas, x_r, x_i, msgs), val);
}

TEST(StanMath_integrate_1d, TestThrows) {
  EXPECT_THROW(stan::math::integrate_1d(f1{}, 1.0, 0.0, std::vector<double>(), {}, {}, msgs, 1e-6), std::domain_error);
  EXPECT_THROW(stan::math::integrate_1d(f1{}, 0.0, std::numeric_limits<double>::quiet_NaN(), std::vector<double>(), {}, {}, msgs, 1e-6), std::domain_error);
  EXPECT_THROW(stan::math::integrate_1d(f1{}, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::vector<double>(), {}, {}, msgs, 1e-6), std::domain_error);
}

TEST(StanMath_integrate_1d, test1) {
  // Tricky integral from Boost docs + limit at infinity
  test_integration(f1{}, 0.0, std::numeric_limits<double>::infinity(), {}, {}, {}, 1.772453850905516);
  // Tricky integral from Boost docs
  test_integration(f2{}, 0.0, 1.0, {}, {}, {}, 1.198140234735592);
  // Zero crossing integral + limit at infinity
  test_integration(f3{}, -2.0, std::numeric_limits<double>::infinity(), {}, {}, {}, 7.38905609893065);
  // Easy integrals
  test_integration(f4{}, 0.2, 0.7, { 0.5 }, {}, {}, 1.04235);
  test_integration(f5{}, -0.2, 0.7, { 0.4, 0.4 }, {}, {}, 1.396622);
  test_integration(f4{}, 0.0, 0.0, { 0.5 }, {}, {}, 0.0);
  test_integration(f5{}, 1.0, 1.0, { 0.4, 0.4 }, {}, {}, 0.0);
  // Test x_i
  test_integration(f6{}, -0.2, 2.9, { 6.0, 5.1 }, {}, {4}, 4131.985);
  // Test x_r
  test_integration(f7{}, -0.2, 2.9, {}, { 4.0, 6.0, 5.1 }, {}, 24219.99);
  // Both limits at infinity + test x_r/x_i
  test_integration(f8{}, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), { 5.0 }, { 1.7 }, { 2 }, 3.013171546539377);
  // Both limits at infinity + test x_i
  test_integration(f9{}, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), { 1.3 }, {}, { 4 }, 2.372032924895055);
}
