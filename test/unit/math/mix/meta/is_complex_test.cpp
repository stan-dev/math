#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>
#include <string>

template <bool expected, typename T>
void expect_is_complex() {
  EXPECT_EQ(expected, stan::is_complex<T>::value);
}

template <typename T>
void test_is_complex() {
  expect_is_complex<true, std::complex<T>>();
  expect_is_complex<true, const std::complex<T>>();
  expect_is_complex<true, std::complex<T>&>();
  expect_is_complex<true, const std::complex<T>&>();
  expect_is_complex<false, std::complex<T>*>();
  expect_is_complex<false, const std::complex<T>*>();
}

TEST(stanMathMix, isComplex) {
  using stan::math::fvar;
  using stan::math::var;
  test_is_complex<double>();
  test_is_complex<var>();
  test_is_complex<fvar<double>>();
  test_is_complex<fvar<fvar<double>>>();
  test_is_complex<fvar<var>>();
  test_is_complex<fvar<fvar<var>>>();

  expect_is_complex<false, bool>();
  expect_is_complex<false, int>();
  expect_is_complex<false, size_t>();
  expect_is_complex<false, double>();
  expect_is_complex<false, std::string>();
}
