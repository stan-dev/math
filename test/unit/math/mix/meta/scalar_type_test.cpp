#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <complex>
#include <vector>

template <typename E, typename T>
void expect_scalar_type() {
  using stan::scalar_type;
  using stan::scalar_type_t;
  using test::expect_same_type;
  expect_same_type<E, typename scalar_type<T>::type>();
  expect_same_type<E, scalar_type_t<T>>();
}

template <typename T>
void test_scalar_type() {
  using Eigen::Matrix;
  using std::complex;
  using std::vector;
  using c_t = complex<T>;

  expect_scalar_type<T, T>();
  expect_scalar_type<T, Matrix<T, -1, -1>>();
  expect_scalar_type<T, Matrix<T, -1, 1>>();
  expect_scalar_type<T, Matrix<T, 1, -1>>();
  expect_scalar_type<T, vector<T>>();

  expect_scalar_type<T&, T&>();
  expect_scalar_type<T, Matrix<T, -1, -1>&>();
  expect_scalar_type<T, Matrix<T, -1, 1>&>();
  expect_scalar_type<T, Matrix<T, 1, -1>&>();
  expect_scalar_type<T, vector<T>&>();

  expect_scalar_type<T, const T>();
  expect_scalar_type<T, const Matrix<T, -1, -1>>();
  expect_scalar_type<T, const Matrix<T, -1, 1>>();
  expect_scalar_type<T, const Matrix<T, 1, -1>>();
  expect_scalar_type<T, const vector<T>>();

  expect_scalar_type<c_t, c_t>();
  expect_scalar_type<c_t, Matrix<c_t, -1, -1>>();
  expect_scalar_type<c_t, Matrix<c_t, -1, 1>>();
  expect_scalar_type<c_t, Matrix<c_t, 1, -1>>();
  expect_scalar_type<c_t, vector<c_t>>();

  expect_scalar_type<c_t, const c_t>();
  expect_scalar_type<c_t, const Matrix<c_t, -1, -1>>();
  expect_scalar_type<c_t, const Matrix<c_t, -1, 1>>();
  expect_scalar_type<c_t, const Matrix<c_t, 1, -1>>();
  expect_scalar_type<c_t, const vector<c_t>>();

  expect_scalar_type<T const*, T const*>();
  expect_scalar_type<T const*, const vector<T const*>>();
  expect_scalar_type<T const*, const vector<vector<T const*>>>();

  expect_scalar_type<c_t, c_t&>();
  expect_scalar_type<c_t, Matrix<c_t, -1, -1>&>();
  expect_scalar_type<c_t, Matrix<c_t, -1, 1>&>();
  expect_scalar_type<c_t, Matrix<c_t, 1, -1>&>();
  expect_scalar_type<c_t, vector<c_t>&>();
}

TEST(MathMetaMix, ScalarTypeInt) {
  using std::vector;
  expect_scalar_type<int, int>();
  expect_scalar_type<int&, int&>();
  expect_scalar_type<int, vector<int>>();
}

TEST(mathMetaMix, scalarType) {
  using stan::math::fvar;
  using stan::math::var;
  test_scalar_type<double>();
  test_scalar_type<var>();
  test_scalar_type<fvar<double>>();
  test_scalar_type<fvar<var>>();
  test_scalar_type<fvar<fvar<double>>>();
  test_scalar_type<fvar<fvar<var>>>();
}
