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

template <bool expected, typename T>
void expect_is_vt_complex() {
  EXPECT_EQ(expected, stan::is_vt_complex<T>::value);
}

template <typename T>
void test_is_vt_complex() {
  using complex_std_vec = std::vector<std::complex<T>>;
  using complex_eigen_vec = Eigen::Matrix<std::complex<T>, -1, 1>;
  expect_is_vt_complex<true, complex_std_vec>();
  expect_is_vt_complex<true, const complex_std_vec>();
  expect_is_vt_complex<true, complex_std_vec&>();
  expect_is_vt_complex<true, const complex_std_vec&>();
  expect_is_vt_complex<false, complex_std_vec*>();
  expect_is_vt_complex<false, const complex_std_vec*>();

  expect_is_vt_complex<true, complex_eigen_vec>();
  expect_is_vt_complex<true, const complex_eigen_vec>();
  expect_is_vt_complex<true, complex_eigen_vec&>();
  expect_is_vt_complex<true, const complex_eigen_vec&>();
  expect_is_vt_complex<false, complex_eigen_vec*>();
  expect_is_vt_complex<false, const complex_eigen_vec*>();
}

TEST(stanMathMix, is_vt_Complex) {
  using stan::math::fvar;
  using stan::math::var;
  test_is_vt_complex<double>();
  test_is_vt_complex<var>();
  test_is_vt_complex<fvar<double>>();
  test_is_vt_complex<fvar<fvar<double>>>();
  test_is_vt_complex<fvar<var>>();
  test_is_vt_complex<fvar<fvar<var>>>();

  expect_is_vt_complex<false, std::vector<bool>>();
  expect_is_vt_complex<false, std::vector<int>>();
  expect_is_vt_complex<false, std::vector<size_t>>();
  expect_is_vt_complex<false, std::vector<double>>();
  expect_is_vt_complex<false, std::vector<std::string>>();
}

template <bool expected, typename T>
void expect_is_vt_not_complex() {
  EXPECT_EQ(expected, stan::is_vt_not_complex<T>::value);
}

template <typename T>
void test_is_vt_not_complex() {
  using complex_std_vec = std::vector<std::complex<T>>;
  using complex_eigen_vec = Eigen::Matrix<std::complex<T>, -1, 1>;
  expect_is_vt_not_complex<!true, complex_std_vec>();
  expect_is_vt_not_complex<!true, const complex_std_vec>();
  expect_is_vt_not_complex<!true, complex_std_vec&>();
  expect_is_vt_not_complex<!true, const complex_std_vec&>();
  expect_is_vt_not_complex<!false, complex_std_vec*>();
  expect_is_vt_not_complex<!false, const complex_std_vec*>();

  expect_is_vt_not_complex<!true, complex_eigen_vec>();
  expect_is_vt_not_complex<!true, const complex_eigen_vec>();
  expect_is_vt_not_complex<!true, complex_eigen_vec&>();
  expect_is_vt_not_complex<!true, const complex_eigen_vec&>();
  expect_is_vt_not_complex<!false, complex_eigen_vec*>();
  expect_is_vt_not_complex<!false, const complex_eigen_vec*>();
}

TEST(stanMathMix, is_vt_not_Complex) {
  using stan::math::fvar;
  using stan::math::var;
  test_is_vt_not_complex<double>();
  test_is_vt_not_complex<var>();
  test_is_vt_not_complex<fvar<double>>();
  test_is_vt_not_complex<fvar<fvar<double>>>();
  test_is_vt_not_complex<fvar<var>>();
  test_is_vt_not_complex<fvar<fvar<var>>>();

  expect_is_vt_not_complex<!false, std::vector<bool>>();
  expect_is_vt_not_complex<!false, std::vector<int>>();
  expect_is_vt_not_complex<!false, std::vector<size_t>>();
  expect_is_vt_not_complex<!false, std::vector<double>>();
  expect_is_vt_not_complex<!false, std::vector<std::string>>();
}
