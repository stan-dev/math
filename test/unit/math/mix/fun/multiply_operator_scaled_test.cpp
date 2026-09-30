#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Scaled matrix products `s * (A * B)`, `(A * B) * s` and `(s * A) * B` are
// evaluated by Stan's product hooks with the scalar folded into the product.

using stan::test::multiply_operator::expect_hooked;
using stan::test::multiply_operator::test_matrix;

TEST(mathMixMultiplyOperator, scaledProductsHooked) {
  const Eigen::MatrixXd b = test_matrix(2, 2, -0.3);
  const Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  double s = 0.7;
  // autodiff scale and left operand, double right operand
  expect_hooked([&b](const auto& s, const auto& a) { return s * (a * b); }, s,
                a);
  expect_hooked([&b](const auto& s, const auto& a) { return (a * b) * s; }, s,
                a);
  expect_hooked([&b](const auto& s, const auto& a) { return (s * a) * b; }, s,
                a);
  // autodiff scale and right operand, double left operand
  expect_hooked([&a](const auto& s, const auto& b) { return s * (a * b); }, s,
                b);
  expect_hooked([&a](const auto& s, const auto& b) { return a * (b * s); }, s,
                b);
}

TEST(mathMixMultiplyOperator, scaledProducts) {
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(2, 2, -0.3);
  double s = 0.7;
  auto f_left = [](const auto& s, const auto& a, const auto& b) {
    return (s * (a * b)).eval();
  };
  auto f_right = [](const auto& s, const auto& a, const auto& b) {
    return ((a * b) * s).eval();
  };
  auto f_inner = [](const auto& s, const auto& a, const auto& b) {
    return ((s * a) * (b * s)).eval();
  };
  auto f_sum = [](const auto& s, const auto& a, const auto& b) {
    return (a.block(0, 0, 2, 2) + s * (a.block(0, 0, 2, 2) * b)).eval();
  };
  stan::test::expect_ad(f_left, s, a, b);
  stan::test::expect_ad(f_right, s, a, b);
  stan::test::expect_ad(f_inner, s, a, b);
  stan::test::expect_ad(f_sum, s, a, b);
}
