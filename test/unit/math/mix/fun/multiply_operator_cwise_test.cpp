#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Values and derivatives of coefficient-wise products and products with a
// scalar of autodiff matrices and arrays. Which Stan function computes each
// form is tested in test/unit/math/rev/core/eigen_product_hook_test.cpp.

using stan::test::multiply_operator::test_matrix;

TEST(mathMixMultiplyOperator, arrayTimesArray) {
  auto f = [](const auto& a, const auto& b) {
    return (a.array() * b.array()).matrix().eval();
  };
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.4);
  stan::test::expect_ad(f, a, b);
}

TEST(mathMixMultiplyOperator, cwiseProduct) {
  auto f
      = [](const auto& a, const auto& b) { return a.cwiseProduct(b).eval(); };
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.4);
  stan::test::expect_ad(f, a, b);
}

TEST(mathMixMultiplyOperator, cwiseNestedAndPartialRead) {
  auto f_nested = [](const auto& a, const auto& b) {
    return ((a.array() * b.array()).exp() * a.array()).matrix().eval();
  };
  auto f_head = [](const auto& a, const auto& b) {
    return (a.array() * b.array()).matrix().col(0).head(2).eval();
  };
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.4);
  stan::test::expect_ad(f_nested, a, b);
  stan::test::expect_ad(f_head, a, b);
}

TEST(mathMixMultiplyOperator, scalarTimesMatrix) {
  auto left = [](const auto& s, const auto& a) { return (s * a).eval(); };
  auto right = [](const auto& s, const auto& a) { return (a * s).eval(); };
  auto array = [](const auto& s, const auto& a) {
    return (a.array() * s).matrix().eval();
  };
  double s = 1.3;
  Eigen::MatrixXd a = test_matrix(2, 3, 0.2);
  stan::test::expect_ad(left, s, a);
  stan::test::expect_ad(right, s, a);
  stan::test::expect_ad(array, s, a);
}

TEST(mathMixMultiplyOperator, matrixTimesInt) {
  auto f = [](const auto& a) { return (a * 3).eval(); };
  stan::test::expect_ad(f, test_matrix(2, 3, 0.2));
}
