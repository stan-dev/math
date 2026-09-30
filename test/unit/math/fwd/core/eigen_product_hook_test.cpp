#include <stan/math/fwd.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>

// Eigen products of fvar<double> matrices go through Stan's product hooks:
// matrix products are split into double products of values and tangents,
// coefficient-wise products stay lazy.

namespace {
using stan::math::fvar;
using fvar_d = fvar<double>;
using stan::math::internal::is_stan_ad_product_hook_v;

template <int R = -1, int C = -1>
Eigen::Matrix<fvar_d, R, C> fvar_mat(int rows, int cols, double offset) {
  Eigen::Matrix<fvar_d, R, C> x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = fvar_d(0.1 * (i + 1) - 0.2 * j + offset, 0.3 * j - 0.1 * i);
    }
  }
  return x;
}

template <typename A, typename B>
void expect_fvar_product(const A& a, const B& b) {
  Eigen::Matrix<fvar_d, -1, -1> c = a * b;
  ASSERT_EQ(c.rows(), a.rows());
  ASSERT_EQ(c.cols(), b.cols());
  for (int j = 0; j < c.cols(); ++j) {
    for (int i = 0; i < c.rows(); ++i) {
      fvar_d expected = 0;
      for (int k = 0; k < a.cols(); ++k) {
        expected += a(i, k) * b(k, j);
      }
      EXPECT_NEAR(c(i, j).val_, expected.val_, 1e-12);
      EXPECT_NEAR(c(i, j).d_, expected.d_, 1e-12);
    }
  }
}
}  // namespace

TEST(mathFwdEigenProductHook, selected) {
  using m = Eigen::Matrix<fvar_d, -1, -1>;
  using md = Eigen::MatrixXd;
  EXPECT_TRUE((is_stan_ad_product_hook_v<std::decay_t<decltype(
                   std::declval<const m&>() * std::declval<const md&>())>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<std::decay_t<decltype(
                   std::declval<const m&>() * std::declval<const m&>())>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<
               std::decay_t<decltype(std::declval<const m&>().array()
                                     * std::declval<const m&>().array())>>));
  EXPECT_FALSE((is_stan_ad_product_hook_v<std::decay_t<decltype(
                    std::declval<const md&>() * std::declval<const md&>())>>));
}

TEST(mathFwdEigenProductHook, matmulValuesAndTangents) {
  expect_fvar_product(fvar_mat(3, 4, 0.1), fvar_mat(4, 2, -0.2));
  const Eigen::MatrixXd d = Eigen::MatrixXd::Constant(4, 3, 0.5);
  expect_fvar_product(fvar_mat(2, 4, 0.1), d);
  expect_fvar_product(d.transpose(), fvar_mat(4, 1, 0.3));
  expect_fvar_product(fvar_mat<2, 2>(2, 2, 0.1), fvar_mat<2, 2>(2, 2, 0.4));
}

TEST(mathFwdEigenProductHook, cwiseIsLazy) {
  auto a = fvar_mat(3, 2, 0.1);
  auto b = fvar_mat(3, 2, -0.4);
  auto xpr = a.array() * b.array();
  Eigen::internal::evaluator<decltype(xpr)> eval(xpr);
  // lazy evaluation: the evaluator does not own a plain result
  EXPECT_FALSE(
      static_cast<bool>(Eigen::internal::evaluator<decltype(xpr)>::Flags
                        & Eigen::EvalBeforeNestingBit));
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 3; ++i) {
      const fvar_d expected = a(i, j) * b(i, j);
      EXPECT_NEAR(eval.coeff(i, j).val_, expected.val_, 1e-12);
      EXPECT_NEAR(eval.coeff(i, j).d_, expected.d_, 1e-12);
    }
  }
}
