#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>
#include <vector>

template <typename R, typename T>
void expect_base() {
  EXPECT_SAME_TYPE(R, typename stan::base_type<T>::type);
  EXPECT_SAME_TYPE(R, stan::base_type_t<T>);
  EXPECT_SAME_TYPE(R, typename stan::base_type<T&>::type);
  EXPECT_SAME_TYPE(R, stan::base_type_t<T&>);
  EXPECT_SAME_TYPE(R, typename stan::base_type<const T&>::type);
  EXPECT_SAME_TYPE(R, stan::base_type_t<const T&>);
  EXPECT_SAME_TYPE(R, typename stan::base_type<const T>::type);
  EXPECT_SAME_TYPE(R, stan::base_type_t<const T>);
}

template <typename T>
void test_base() {
  // scalar types
  expect_base<T, T>();

  // array types
  expect_base<T, std::vector<T>>();
  expect_base<T, std::vector<std::vector<T>>>();

  // matrix types
  expect_base<T, Eigen::Matrix<T, -1, -1>>();
  expect_base<T, Eigen::Matrix<T, -1, 1>>();
  expect_base<T, Eigen::Matrix<T, 1, -1>>();

  // complex types
  expect_base<T, std::complex<T>>();

  // higher-level containers
  expect_base<T, std::vector<Eigen::Matrix<T, -1, -1>>>();
  expect_base<T, std::vector<Eigen::Matrix<T, -1, 1>>>();
  expect_base<T, std::vector<Eigen::Matrix<T, 1, -1>>>();
  expect_base<T, std::vector<std::complex<T>>>();
}

TEST(mathMetaMix, baseType) {
  using stan::math::fvar;
  using stan::math::var;
  // 1-arg special cases where result is min double
  expect_base<float, float>();
  expect_base<int, int>();
  expect_base<int, std::vector<int>>();
  expect_base<int, std::vector<std::vector<int>>>();

  // cases where result is given real type
  test_base<float>();
  test_base<double>();
  test_base<var>();
  test_base<fvar<double>>();
  test_base<fvar<fvar<double>>>();
  test_base<fvar<var>>();
  test_base<fvar<fvar<var>>>();
}
