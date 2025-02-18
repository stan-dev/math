#include <test/unit/math/test_ad.hpp>
#include <limits>

/**
 * NOTE: We do not test values of (0.0, 0.0) as inputs.
 *  This is because the testing framework uses
 *  finite difference as the comparison. The small finite
 *  purtubations lead to cases where the value of the
 *  function inputs is either slightly negative or slightly positive,
 *  which can lead to the function returning NaN when we would expect a
 *  value of 0.
 */
TEST(mathMixScalFun, multiplyLog2_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::multiply_log;
    return multiply_log(x1, x2);
  };
  auto expect_test = [](auto&& f, auto&& x1, auto&& x2) {
    stan::test::expect_ad(f, x1, x2);
    stan::test::expect_ad(f, x2, x1);
    stan::test::expect_ad_matvar(f, x1, x2);
    stan::test::expect_ad_matvar(f, x2, x1);
  };
  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;

  Eigen::VectorXd x1(3);
  x1 << 1.0, 2.0, 3.0;
  Eigen::RowVectorXd x2(3);
  x2 << 1.0, 2.0, 3.0;
  Eigen::MatrixXd x3(2, 3);
  x3 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  stan::test::expect_ad(f, x1, x1);
  stan::test::expect_ad(f, x2, x2);
  stan::test::expect_ad(f, x3, x3);
  stan::test::expect_ad_matvar(f, x1, x1);
  stan::test::expect_ad_matvar(f, x2, x2);
  stan::test::expect_ad_matvar(f, x3, x3);
  expect_test(f, x1, 2.0);
  expect_test(f, x2, 2.5);
  expect_test(f, x3, 5.5);

  Eigen::VectorXd x4(0);
  Eigen::RowVectorXd x5(0);
  Eigen::MatrixXd x6(0, 0);

  stan::test::expect_ad(f, x4, x4);
  stan::test::expect_ad(f, x5, x5);
  stan::test::expect_ad(f, x6, x6);
  stan::test::expect_ad_matvar(f, x4, x4);
  stan::test::expect_ad_matvar(f, x5, x5);
  stan::test::expect_ad_matvar(f, x6, x6);
  expect_test(f, x4, 2.0);
  expect_test(f, x5, 3.1);
  expect_test(f, x6, 5.5);
}

TEST(mathMixScalFun, multiplyLog2_zero_vec_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::multiply_log;
    return multiply_log(x1, x2);
  };
  auto expect_test = [](auto&& f, auto&& x1, auto&& x2) {
    stan::test::expect_ad(f, x1, x2);
    stan::test::expect_ad(f, x2, x1);
    stan::test::expect_ad_matvar(f, x1, x2);
    stan::test::expect_ad_matvar(f, x2, x1);
  };

  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd x1(3);
  x1 << 1.0, 2.0, 3.0;
  expect_test(f, zero_vec, x1);
}
TEST(mathMixScalFun, multiplyLog2_zero_vec_scalar) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::multiply_log;
    return multiply_log(x1, x2);
  };

  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd x1(3);
  x1 << 1.0, 2.0, 3.0;
  auto expect_test = [](auto&& f, auto&& x1, auto&& x2) {
    stan::test::expect_ad(f, x1, x2);
    stan::test::expect_ad(f, x2, x1);
    stan::test::expect_ad_matvar(f, x1, x2);
    stan::test::expect_ad_matvar(f, x2, x1);
  };
  expect_test(f, x1, 0.0);
}
