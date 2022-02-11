#include <test/unit/math/test_ad.hpp>
#include <vector>
#include <gtest/gtest.h>

TEST(mathMixMatFun, gen_inv_test_1) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd t(0, 0);
  expect_ad(f, t);
  expect_ad_matvar(f, t);
}

TEST(mathMixMatFun, gen_inv_test_2) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd u(1, 1);
  u << 2;
  expect_ad(f, u);
  expect_ad_matvar(f, u);
}

TEST(mathMixMatFun, gen_inv_test_3) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd v(2, 3);
  v << 1, 3, 5, 2, 4, 6;
  expect_ad(f, v);
  expect_ad_matvar(f, v);
}

TEST(mathMixMatFun, gen_inv_test_4) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd v(2, 3);
  v << 1.9, 1.3, 2.5, 0.4, 1.7, 0.1;
  expect_ad(f, v);
  expect_ad_matvar(f, v);
}

TEST(mathMixMatFun, gen_inv_test_5) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd s(2, 4);
  s << 3.4, 2, 5, 1.2, 2, 1, 3.2, 3.1;
  expect_ad(f, s);
  expect_ad_matvar(f, s);
}

TEST(mathMixMatFun, gen_inv_test_6) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  // issues around zero require looser tolerances for hessians
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 2.0;
  tols.hessian_fvar_hessian_ = 2.0;

  Eigen::MatrixXd w(3, 4);
  w << 2, 3, 5, 7, 11, 13, 17, 19, 23, 25, 27, 29;
  expect_ad(tols, f, w);
  expect_ad_matvar(f, w);
}

TEST(mathMixMatFun, gen_inv_test_8) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  Eigen::MatrixXd z(2, 2);
  z << 1, 2, 5, std::numeric_limits<double>::quiet_NaN();
  EXPECT_NO_THROW(stan::math::generalized_inverse(z));
}

TEST(mathMixMatFun, gen_inv_test_9) {
  using stan::test::expect_ad;
  using stan::test::expect_ad_matvar;

  auto f = [](const auto& G) { return stan::math::generalized_inverse(G); };

  // autodiff throws, so following fails (throw behavior must match to pass)

  Eigen::MatrixXd a(2, 2);
  a << 1.9, 0.3, 0.3, std::numeric_limits<double>::infinity();
  expect_ad(f, a);
  expect_ad_matvar(f, a);
}
