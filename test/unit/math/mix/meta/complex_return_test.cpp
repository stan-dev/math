#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>

template <typename T>
void test_complex_return() {
  using stan::return_type;
  using stan::math::complex_return_t;
  using test::expect_same_type;
  using c_t = std::complex<T>;
  using c_d = std::complex<double>;

  // unary
  expect_same_type<c_t, complex_return_t<c_t>>();

  // binary
  expect_same_type<c_t, complex_return_t<c_t, double>>();
  expect_same_type<c_t, complex_return_t<double, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, int>>();
  expect_same_type<c_t, complex_return_t<int, c_t>>();

  expect_same_type<c_t, complex_return_t<c_t, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, int>>();
  expect_same_type<c_t, complex_return_t<int, c_t>>();

  // ternary
  expect_same_type<c_t, complex_return_t<c_t, double, double>>();
  expect_same_type<c_t, complex_return_t<double, c_t, double>>();
  expect_same_type<c_t, complex_return_t<double, double, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, c_t, double>>();
  expect_same_type<c_t, complex_return_t<double, c_t, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, double, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, c_t, c_t>>();

  expect_same_type<c_t, complex_return_t<c_t, c_d, double>>();
  expect_same_type<c_t, complex_return_t<c_t, double, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, c_t, double>>();
  expect_same_type<c_t, complex_return_t<double, c_t, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, double, c_t>>();
  expect_same_type<c_t, complex_return_t<double, c_d, c_t>>();

  expect_same_type<c_t, complex_return_t<c_t, c_d, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, c_t, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, c_d, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, c_t, c_d>>();
  expect_same_type<c_t, complex_return_t<c_d, c_t, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, c_d, c_t>>();
  expect_same_type<c_t, complex_return_t<c_t, c_t, c_t>>();
}
TEST(mathMixFun, valueOfRecComplex) {
  using stan::math::fvar;
  using stan::math::var;
  test_complex_return<double>();
  test_complex_return<var>();
  test_complex_return<fvar<double>>();
  test_complex_return<fvar<fvar<double>>>();
  test_complex_return<fvar<var>>();
  test_complex_return<fvar<fvar<var>>>();
}
