#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(test_unit_math_test_ad, test_ad_unary) {
  Eigen::MatrixXd x(2, 2);
  x << 1.9, 0.3, 0.3, 1.7;

  auto g = [](const auto& u) { return stan::math::inverse(u); };
  stan::test::expect_ad(g, x);

  // Note:  the above is how tests need to be written;  the
  // functor is required because the following won't compile
  // because the template parameter F for expect_ad can't be deduced
  //   stan::test::expect_ad(inverse, x);
}

TEST(test_unit_math_test_ad, test_ad_binary) {
  auto g
      = [](const auto& u, const auto& v) { return stan::math::multiply(u, v); };

  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  double y = -3;
  stan::test::expect_ad(g, x, y);

  double y2 = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(g, x, y2);
}

// VECTORIZED UNARY FUNCTION THAT PASSES

// log10 is vectorized, so uses vectorized test
TEST(test_unit_math_test_ad, expect_ad_vectorized) {
  auto g = [](const auto& u) { return stan::math::log10(u); };

  stan::test::expect_ad_vectorized(g, 3.2);
}

// OVERLOAD THAT PASSES

// double overload matches template behavior and passes tests
template <typename T>
T f_match(const T& x) {
  return -2 * x;
}
double f_match(const double& x) { return -2 * x; }
TEST(test_ad, match) {
  double x = 3.2;
  auto g = [](const auto& u) { return f_match(u); };
  stan::test::expect_ad(g, x);
}

// OVERLOAD THAT FAILS DUE TO MISMATCHED VALUE / DERIVATIVE

// double overload does not match template behavior for value
template <typename T>
T f_mismatch(const T& x) {
  return -2 * x;
}
stan::math::var f_mismatch(const stan::math::var& x) { return 2 * x; }
TEST(test_ad, mismatch) {
  double x = 3.2;
  auto g = [](const auto& u) { return f_mismatch(u); };
  // include following line to show exception error behavior
  // stan::test::expect_ad(g, x);
}

// OVERLOAD THAT FAILS DUE TO MISMATCHED EXCEPTION CONDITIONS

// double overload does not match template behavior for exceptionsalue
template <typename T>
T f_misthrow(const T& x) {
  return -2 * x;
}
double f_misthrow(const double& x) {
  throw std::runtime_error("f_misthrow(double) called");
  return -2 * x;
}

TEST(test_ad, misthrow) {
  double x = 1.73;
  auto h = [](const auto& u) { return f_misthrow(u); };
  // include following line to show exception error behavior
  // stan::test::expect_ad(h, x);
}

template <typename T>
T foo(const T& x) {
  // std::cout << "called foo(T)" << std::endl;
  return x / 2;
}

double foo(double x) {
  // std::cout << "called foo(double)" << std::endl;
  return x / 2;
}

double foo(int x) {
  // std::cout << "called foo(int)" << std::endl;
  // return x / 2.0; // passes
  // return x / 2;   // fails on odd numbers
  return 2 * x;  // always fails
}

TEST(test_ad, integerGetsPassed) {
  auto h = [](const auto& u) { return foo(u); };
  stan::test::expect_ad(h, 3.0);  // passes for double
  // stan::test::expect_ad(h, 3);    // fails
}
