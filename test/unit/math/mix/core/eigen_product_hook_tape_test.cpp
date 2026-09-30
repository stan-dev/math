#include <stan/math/mix.hpp>
#include <test/unit/math/ad_tape_util.hpp>
#include <gtest/gtest.h>
#include <type_traits>

// Reverse-mode nodes recorded by Eigen products of fvar<var> matrices,
// which Stan splits into products of their var values and tangents.

namespace {
using stan::math::fvar;
using stan::math::var;
using fvar_v = fvar<var>;

Eigen::Matrix<fvar_v, -1, -1> fvar_var_mat(int rows, int cols, double offset) {
  Eigen::Matrix<fvar_v, -1, -1> x(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      x(i, j) = fvar_v(var(0.1 * (i + 1) - 0.2 * j + offset), var(1.0));
    }
  }
  return x;
}
}  // namespace

TEST(mathMixEigenProductHook, fvarVarProductsAreHooked) {
  using stan::math::internal::is_stan_ad_product_hook_v;
  using m = Eigen::Matrix<fvar_v, -1, -1>;
  using md = Eigen::MatrixXd;
  EXPECT_TRUE((is_stan_ad_product_hook_v<std::decay_t<decltype(
                   std::declval<const m&>() * std::declval<const md&>())>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<std::decay_t<decltype(
                   std::declval<const m&>() * std::declval<const m&>())>>));
  EXPECT_TRUE((is_stan_ad_product_hook_v<
               std::decay_t<decltype(std::declval<const m&>().array()
                                     * std::declval<const m&>().array())>>));
}

TEST(mathMixEigenProductHook, fvarVarMatmulTape) {
  stan::math::nested_rev_autodiff nested;
  const int n = 6, k = 5, m = 4;
  auto a = fvar_var_mat(n, k, 0.1);
  auto b = fvar_var_mat(k, m, -0.3);
  const stan::test::tape_delta delta;
  Eigen::Matrix<fvar_v, -1, -1> c = a * b;
  // one callback for the values and one for the concatenated tangent product
  EXPECT_EQ(delta.chain(), 2);
  EXPECT_EQ(delta.nochain(), 2 * n * m);
}

TEST(mathMixEigenProductHook, fvarVarDoubleMatmulTape) {
  stan::math::nested_rev_autodiff nested;
  const int n = 6, k = 5;
  const Eigen::MatrixXd x = Eigen::MatrixXd::Constant(n, k, 0.25);
  auto beta = fvar_var_mat(k, 1, 0.2);
  const stan::test::tape_delta delta;
  Eigen::Matrix<fvar_v, -1, 1> xb = x * beta;
  EXPECT_EQ(delta.chain(), 2);
  EXPECT_EQ(delta.nochain(), 2 * n);
}

TEST(mathMixEigenProductHook, fvarVarCwiseTape) {
  stan::math::nested_rev_autodiff nested;
  const int n = 7;
  auto a = fvar_var_mat(n, 1, 0.1);
  auto b = fvar_var_mat(n, 1, -0.3);
  const stan::test::tape_delta delta;
  Eigen::Matrix<fvar_v, -1, 1> c = (a.array() * b.array()).matrix();
  // values: one callback; tangent: two callbacks plus one vari per
  // coefficient for their sum
  EXPECT_EQ(delta.chain(), 3 + n);
  EXPECT_EQ(delta.nochain(), 3 * n);
}
