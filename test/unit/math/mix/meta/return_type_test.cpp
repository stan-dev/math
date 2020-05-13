#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>
#include <vector>

using stan::return_type;
using stan::math::fvar;
using stan::math::var;
using d_t = double;
using v_t = var;
using fd_t = fvar<double>;
using ffd_t = fvar<fd_t>;
using fv_t = fvar<var>;
using ffv_t = fvar<fv_t>;
using md_t = Eigen::MatrixXd;
using vd_t = Eigen::VectorXd;
using rvd_t = Eigen::RowVectorXd;

template <typename R, typename... Ts>
void expect_return() {
  test::expect_same_type<R, typename stan::return_type<Ts...>::type>();
  test::expect_same_type<R, stan::return_type_t<Ts...>>();
}

template <typename T>
void test_return() {
  // scalar types
  expect_return<T, T>();
  expect_return<T, T, int>();
  expect_return<T, int, T>();
  expect_return<T, d_t, T, d_t, int, d_t, float, float, float, T, int>();

  // array types
  expect_return<T, std::vector<T>>();
  expect_return<T, std::vector<T>, int>();
  expect_return<T, double, std::vector<T>>();

  // matrix types
  expect_return<T, Eigen::Matrix<T, -1, -1>>();
  expect_return<T, Eigen::Matrix<T, -1, 1>>();
  expect_return<T, Eigen::Matrix<T, -1, 1>, std::vector<double>>();
  expect_return<T, Eigen::Matrix<T, 1, -1>>();
  expect_return<T, Eigen::Matrix<T, 1, -1>, T>();
  expect_return<T, T, Eigen::Matrix<T, 1, -1>, int>();
  expect_return<T, double, Eigen::Matrix<T, 1, -1>>();
  expect_return<T, Eigen::Matrix<T, 1, -1>, double>();
  expect_return<T, Eigen::Matrix<T, 1, -1>, int, Eigen::Matrix<T, -1, -1>>();
  expect_return<T, Eigen::Matrix<T, 1, -1>, int,
                std::vector<Eigen::Matrix<double, -1, -1>>>();

  // complex types
  expect_return<std::complex<T>, std::complex<T>>();

  expect_return<std::complex<T>, int, std::complex<T>>();
  expect_return<std::complex<T>, std::complex<T>, int>();

  expect_return<std::complex<T>, double, std::complex<T>>();
  expect_return<std::complex<T>, std::complex<T>, double>();

  expect_return<std::complex<T>, std::complex<double>, std::complex<T>>();
  expect_return<std::complex<T>, std::complex<T>, std::complex<double>>();

  expect_return<std::complex<T>, T, std::complex<T>>();
  expect_return<std::complex<T>, std::complex<T>, T>();

  expect_return<std::complex<T>, std::complex<T>, std::complex<T>>();
  expect_return<std::complex<T>, std::complex<T>, std::complex<T>, T>();
}

TEST(mathMetaMix, returnType) {
  // no-arg case
  expect_return<double>();

  // 1-arg special cases where result is min double
  expect_return<double, float>();
  expect_return<double, int>();
  expect_return<double, std::vector<int>>();
  expect_return<double, Eigen::Matrix<float, -1, 1>>();

  // cases where result is given real type
  test_return<d_t>();
  test_return<v_t>();
  test_return<fd_t>();
  test_return<ffd_t>();
  test_return<fv_t>();
  test_return<ffv_t>();
}
