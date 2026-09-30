#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Values and derivatives of scaled matrix products `s * (A * B)`,
// `(A * B) * s` and `(s * A) * B`. Which Stan function computes each form is
// tested in test/unit/math/rev/core/eigen_product_hook_test.cpp.

using stan::test::multiply_operator::test_matrix;

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
  // Eigen cannot multiply fvar and double matrices until fvar products are
  // routed to Stan's functions, so this test runs in reverse mode only.
  constexpr bool reverse_only = true;
  stan::test::expect_ad<reverse_only>(f_left, s, a, b);
  stan::test::expect_ad<reverse_only>(f_right, s, a, b);
  stan::test::expect_ad<reverse_only>(f_inner, s, a, b);
  stan::test::expect_ad<reverse_only>(f_sum, s, a, b);
}
