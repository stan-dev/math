#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

template <typename R, typename T>
void expect_scalar_type() {
  using stan::scalar_type_t;
  static_assert(std::is_same<R, scalar_type_t<T>>::value, "NOT SAME");
  static_assert(std::is_same<R, scalar_type_t<T&>>::value, "NOT SAME");
  EXPECT_SAME_TYPE(R, scalar_type_t<T>);
  EXPECT_SAME_TYPE(R, scalar_type_t<T>);
  EXPECT_SAME_TYPE(R, scalar_type_t<T&>);
  EXPECT_SAME_TYPE(R, scalar_type_t<T&>);
  EXPECT_SAME_TYPE(R, scalar_type_t<const T&>);
  EXPECT_SAME_TYPE(R, scalar_type_t<const T&>);
  EXPECT_SAME_TYPE(R, scalar_type_t<const T>);
  EXPECT_SAME_TYPE(R, scalar_type_t<const T>);
}

template <typename T>
void test_scalar_type() {
  using Eigen::Matrix;
  using std::complex;
  using std::vector;
  expect_scalar_type<T, T>();
  expect_scalar_type<T, vector<T>>();
  expect_scalar_type<T, vector<vector<T>>>();

  expect_scalar_type<T, Matrix<T, -1, -1>>();
  expect_scalar_type<T, vector<Matrix<T, -1, -1>>>();

  expect_scalar_type<T, Matrix<T, 1, -1>>();
  expect_scalar_type<T, vector<Matrix<T, 1, -1>>>();

  expect_scalar_type<T, Matrix<T, -1, 1>>();
  expect_scalar_type<T, vector<Matrix<T, -1, 1>>>();

  expect_scalar_type<complex<T>, complex<T>>();
  expect_scalar_type<complex<T>, vector<complex<T>>>();

  expect_scalar_type<complex<T>, Matrix<complex<T>, -1, -1>>();
}

TEST(MathMetaPrim, scalarType) {
  using stan::base_type;
  using stan::math::fvar;
  using stan::math::var;
  using d_t = double;
  using v_t = var;
  using fd_t = fvar<double>;
  using ffd_t = fvar<fd_t>;
  using fv_t = fvar<var>;
  using ffv_t = fvar<fv_t>;
  EXPECT_SAME_TYPE(int, stan::scalar_type<int>::type);
  EXPECT_SAME_TYPE(int, stan::scalar_type<std::vector<int>>::type);

  test_scalar_type<d_t>();
  test_scalar_type<d_t const*>();
  test_scalar_type<v_t>();
  test_scalar_type<v_t const*>();
  test_scalar_type<fd_t>();
  test_scalar_type<ffd_t>();
  test_scalar_type<fv_t>();
  test_scalar_type<ffv_t>();
}

TEST(MathMetaPrim, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;
  // These ones work because pointers are always copied
  EXPECT_SAME_TYPE(double const*,
                   scalar_type<const vector<double const*>>::type);
  EXPECT_SAME_TYPE(int const*, scalar_type<const vector<int const*>>::type);
  EXPECT_SAME_TYPE(double const*,
                   scalar_type<const vector<vector<double const*>>>::type);
}
