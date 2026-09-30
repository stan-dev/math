#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Coefficient-wise products and products with a scalar of autodiff
// matrices and arrays are evaluated by Stan's product hooks.

using stan::test::multiply_operator::expect_hooked;
using stan::test::multiply_operator::test_matrix;

TEST(mathMixMultiplyOperator, arrayTimesArray) {
  auto xpr = [](const auto& a, const auto& b) { return a.array() * b.array(); };
  auto f = [](const auto& a, const auto& b) {
    return (a.array() * b.array()).matrix().eval();
  };
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.4);
  expect_hooked(xpr, a, b);
  stan::test::expect_ad(f, a, b);
}

TEST(mathMixMultiplyOperator, cwiseProduct) {
  auto xpr = [](const auto& a, const auto& b) { return a.cwiseProduct(b); };
  auto f
      = [](const auto& a, const auto& b) { return a.cwiseProduct(b).eval(); };
  Eigen::MatrixXd a = test_matrix(3, 2, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.4);
  expect_hooked(xpr, a, b);
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
  auto left = [](const auto& s, const auto& a) { return s * a; };
  auto right = [](const auto& s, const auto& a) { return a * s; };
  auto array = [](const auto& s, const auto& a) { return a.array() * s; };
  double s = 1.3;
  Eigen::MatrixXd a = test_matrix(2, 3, 0.2);
  expect_hooked(left, s, a);
  expect_hooked(right, s, a);
  expect_hooked(array, s, a);
  auto eval = [](const auto& f) {
    return [f](const auto& s, const auto& a) {
      return Eigen::Matrix<stan::return_type_t<decltype(s), decltype(a)>, -1,
                           -1>(f(s, a));
    };
  };
  stan::test::expect_ad(eval(left), s, a);
  stan::test::expect_ad(eval(right), s, a);
  stan::test::expect_ad(eval(array), s, a);
}

TEST(mathMixMultiplyOperator, matrixTimesInt) {
  Eigen::MatrixXd a = test_matrix(2, 3, 0.2);
  // only the matrix can be autodiff here, so check the left promotion
  using stan::math::fvar;
  using stan::math::var;
  EXPECT_TRUE(
      (stan::math::internal::is_stan_ad_product_hook_v<std::decay_t<decltype(
           std::declval<const Eigen::Matrix<var, -1, -1>&>() * 3)>>));
  EXPECT_TRUE(
      (stan::math::internal::is_stan_ad_product_hook_v<std::decay_t<decltype(
           std::declval<const Eigen::Matrix<fvar<var>, -1, -1>&>() * 3)>>));
  auto f = [](const auto& a) { return (a * 3).eval(); };
  stan::test::expect_ad(f, a);
}
