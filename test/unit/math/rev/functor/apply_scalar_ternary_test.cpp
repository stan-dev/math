#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace apply_scalar_ternary_test {
struct fma_functor {
  template <typename T1, typename T2, typename T3>
  auto operator()(const T1& x, const T2& y, const T3& z) const {
    return x * y + z;
  }
};

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
Eigen::ArrayXXd broadcast(const T& x, Eigen::Index rows, Eigen::Index cols) {
  return Eigen::ArrayXXd::Constant(rows, cols, stan::math::value_of(x));
}

template <typename T, stan::require_eigen_t<T>* = nullptr>
Eigen::ArrayXXd broadcast(const T& x, Eigen::Index rows, Eigen::Index cols) {
  return stan::math::value_of(x).array();
}

template <typename T1, typename T2, typename T3>
void expect_fma(const T1& x, const T2& y, const T3& z) {
  Eigen::MatrixXd res = stan::math::value_of(
      stan::math::apply_scalar_ternary(fma_functor(), x, y, z));
  Eigen::MatrixXd expected = broadcast(x, res.rows(), res.cols())
                                 * broadcast(y, res.rows(), res.cols())
                             + broadcast(z, res.rows(), res.cols());
  EXPECT_MATRIX_FLOAT_EQ(expected, res);
}

template <typename T>
void test_arena(const Eigen::MatrixXd& x_val, const Eigen::MatrixXd& y_val,
                const Eigen::MatrixXd& z_val) {
  using stan::math::arena_matrix;
  using scalar_t = typename T::Scalar;
  arena_matrix<T> x = x_val;
  arena_matrix<T> y = y_val;
  arena_matrix<T> z = z_val;
  T z_plain = z_val;
  scalar_t s = 1.5;

  expect_fma(x, y, z);
  expect_fma(x, y, z_plain);
  expect_fma(z_plain, x, y);
  expect_fma(x, y, z_plain.array());
  expect_fma(x, y, s);
  expect_fma(x, s, z);
  expect_fma(s, y, z);
  expect_fma(x, s, s);
  expect_fma(s, y, s);
  expect_fma(s, s, z);

  // rvalue arena_matrix inputs
  Eigen::MatrixXd res = stan::math::value_of(stan::math::apply_scalar_ternary(
      fma_functor(), arena_matrix<T>(x), arena_matrix<T>(y),
      arena_matrix<T>(z)));
  EXPECT_MATRIX_FLOAT_EQ((x_val.array() * y_val.array() + z_val.array()), res);
}
}  // namespace apply_scalar_ternary_test

TEST(AgradRevFunctor, apply_scalar_ternary_arena_matrix) {
  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd y = x.array() - 2.5;
  Eigen::MatrixXd z = x.array() * 0.5;

  apply_scalar_ternary_test::test_arena<Eigen::MatrixXd>(x, y, z);
  apply_scalar_ternary_test::test_arena<stan::math::matrix_v>(x, y, z);
  stan::math::recover_memory();
}

TEST(AgradRevFunctor, apply_scalar_ternary_arena_matrix_std_vector) {
  using stan::math::arena_matrix;
  Eigen::VectorXd x(3);
  x << 1, 2, 3;
  std::vector<arena_matrix<Eigen::VectorXd>> x_vec{x, x};
  std::vector<arena_matrix<stan::math::vector_v>> y_vec{x, x};
  Eigen::VectorXd expected = x.array() * x.array() + 1.5;

  auto res = stan::math::apply_scalar_ternary(
      apply_scalar_ternary_test::fma_functor(), x_vec, y_vec, 1.5);
  auto res_all = stan::math::apply_scalar_ternary(
      apply_scalar_ternary_test::fma_functor(), x_vec, y_vec, x_vec);
  for (size_t i = 0; i < res.size(); ++i) {
    EXPECT_MATRIX_FLOAT_EQ(expected, stan::math::value_of(res[i]));
    EXPECT_MATRIX_FLOAT_EQ((x.array() * x.array() + x.array()),
                           stan::math::value_of(res_all[i]));
  }
  stan::math::recover_memory();
}
