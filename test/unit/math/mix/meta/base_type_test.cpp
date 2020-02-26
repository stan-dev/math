#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>
#include <vector>

using stan::base_type;
using stan::math::fvar;
using stan::math::var;
using d_t = double;
using v_t = var;
using fd_t = fvar<double>;
using ffd_t = fvar<fd_t>;
using fv_t = fvar<var>;
using ffv_t = fvar<fv_t>;

template <typename R, typename T>
void expect_base() {
  test::expect_same_type<R, typename stan::base_type<T>::type>();
  test::expect_same_type<R, stan::base_type_t<T>>();
  test::expect_same_type<R, typename stan::base_type<T&>::type>();
  test::expect_same_type<R, stan::base_type_t<T&>>();
  test::expect_same_type<R, typename stan::base_type<const T&>::type>();
  test::expect_same_type<R, stan::base_type_t<const T&>>();
  test::expect_same_type<R, typename stan::base_type<const T>::type>();
  test::expect_same_type<R, stan::base_type_t<const T>>();
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
  // 1-arg special cases where result is min double
  expect_base<float, float>();
  expect_base<int, int>();
  expect_base<int, std::vector<int>>();
  expect_base<int, std::vector<std::vector<int>>>();

  // cases where result is given real type
  test_base<float>();
  test_base<d_t>();
  test_base<v_t>();
  test_base<fd_t>();
  test_base<ffd_t>();
  test_base<fv_t>();
  test_base<ffv_t>();
}
