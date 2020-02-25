#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
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
  test::expect_same_type<R, typename stan::scalar_type<T>::type>();
  test::expect_same_type<R, stan::scalar_type_t<T>>();
  test::expect_same_type<R, typename stan::scalar_type<T&>::type>();
  test::expect_same_type<R, stan::scalar_type_t<T&>>();
  test::expect_same_type<R, typename stan::scalar_type<const T&>::type>();
  test::expect_same_type<R, stan::scalar_type_t<const T&>>();
  test::expect_same_type<R, typename stan::scalar_type<const T>::type>();
  test::expect_same_type<R, stan::scalar_type_t<const T>>();
}

template <typename T>
void test_base() {
  using Eigen::Matrix;
  using std::complex;
  using std::vector;
  expect_base<T, T>();
  expect_base<T, vector<T>>();
  expect_base<T, vector<const T>>();
  // expect_base<T, vector<T&>>();
  expect_base<T, vector<vector<T>>>();
  expect_base<T, Matrix<T, -1, -1>>();
  expect_base<T, Matrix<T, 1, -1>>();
  expect_base<T, Matrix<T, -1, 1>>();
  expect_base<T, vector<Matrix<T, -1, -1>>>();
  expect_base<T, vector<Matrix<T, 1, -1>>>();
  expect_base<T, vector<Matrix<T, -1, 1>>>();

  expect_base<complex<T>, complex<T>>();
}

TEST(MathMetaPrim, scalarType) {
  test::expect_same_type<int, stan::scalar_type<int>::type>();
  test::expect_same_type<int, stan::scalar_type<std::vector<int>>::type>();

  test_base<d_t>();
  test_base<v_t>();
  test_base<fd_t>();
  test_base<ffd_t>();
  test_base<fv_t>();
  test_base<ffv_t>();
}

TEST(MathMetaPrim, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double const*,
                         scalar_type<const vector<double const*>>::type>();
  test::expect_same_type<int const*,
                         scalar_type<const vector<int const*>>::type>();
  test::expect_same_type<
      double const*, scalar_type<const vector<vector<double const*>>>::type>();
}
