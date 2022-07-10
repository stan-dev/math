#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <gtest/gtest.h>
#include <string>

template <typename T_y, typename T_dof, typename T_scale>
void expect_propto_inv_wishart_cholesky_lpdf(T_y L_Y1, T_dof nu1, T_scale L_S1,
                                             T_y L_Y2, T_dof nu2, T_scale L_S2,
                                             std::string message) {
  expect_eq_diffs(stan::math::inv_wishart_cholesky_lpdf<false>(L_Y1, nu1, L_S1),
                  stan::math::inv_wishart_cholesky_lpdf<false>(L_Y2, nu2, L_S2),
                  stan::math::inv_wishart_cholesky_lpdf<true>(L_Y1, nu1, L_S1),
                  stan::math::inv_wishart_cholesky_lpdf<true>(L_Y2, nu2, L_S2),
                  message);
}

using Eigen::Dynamic;
using Eigen::Matrix;

class AgradDistributionsInvWishartCholesky : public ::testing::Test {
 protected:
  virtual void SetUp() {
    Eigen::MatrixXd Y1(2, 2);
    Eigen::MatrixXd Y2(2, 2);
    Y1 << 2.011108, -11.20661, -11.20661, 112.94139;
    Y2 << 13.4, 12.2, 12.2, 11.5;

    L_Y1 = Y1.llt().matrixL();
    L_Y2 = Y2.llt().matrixL();

    nu1 = 3;
    nu2 = 5.3;

    Eigen::MatrixXd S1(2, 2);
    Eigen::MatrixXd S2(2, 2);
    S1 << 1.848220, 1.899623, 1.899623, 12.751941;
    S2 << 3.0, 1.4, 1.4, 7.0;

    L_S1 = S1.llt().matrixL();
    L_S2 = S2.llt().matrixL();
  }
  Eigen::MatrixXd L_Y1;
  Eigen::MatrixXd L_Y2;
  double nu1;
  double nu2;
  Eigen::MatrixXd L_S1;
  Eigen::MatrixXd L_S2;
};

TEST_F(AgradDistributionsInvWishartCholesky, Propto) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(
      to_var(L_Y1), to_var(nu1), to_var(L_S1), to_var(L_Y2), to_var(nu2),
      to_var(L_S2), "var: L_Y, nu, and L_S");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoY) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(to_var(L_Y1), nu1, L_S1, to_var(L_Y2),
                                          nu1, L_S1, "var: L_Y");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoYNu) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(to_var(L_Y1), to_var(nu1), L_S1,
                                          to_var(L_Y2), to_var(nu2), L_S1,
                                          "var: L_y and nu");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoYLSigma) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(to_var(L_Y1), nu1, to_var(L_S1),
                                          to_var(L_Y2), nu1, to_var(L_S2),
                                          "var: L_Y and L_S");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoNu) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(L_Y1, to_var(nu1), L_S1, L_Y1,
                                          to_var(nu2), L_S1, "var: nu");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoNuL_S) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(L_Y1, to_var(nu1), to_var(L_S1), L_Y1,
                                          to_var(nu2), to_var(L_S2),
                                          "var: nu and sigma");

  stan::math::recover_memory();
}
TEST_F(AgradDistributionsInvWishartCholesky, ProptoL_S) {
  using stan::math::to_var;
  expect_propto_inv_wishart_cholesky_lpdf(L_Y1, nu1, to_var(L_S1), L_Y1, nu1,
                                          to_var(L_S2), "var: L_S");

  stan::math::recover_memory();
}

TEST(InvWishartCholesky, check_varis_on_stack) {
  using stan::math::to_var;
  Eigen::MatrixXd Y(2, 2);
  Y << 2.011108, -11.20661, -11.20661, 112.94139;
  Eigen::MatrixXd L_Y = Y.llt().matrixL();
  double nu = 3;

  Eigen::MatrixXd S(2, 2);
  S << 1.848220, 1.899623, 1.899623, 12.751941;
  Eigen::MatrixXd L_S = S.llt().matrixL();

  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<false>(
      to_var(L_Y), to_var(nu), to_var(L_S)));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<false>(
      to_var(L_Y), to_var(nu), L_S));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<false>(
      to_var(L_Y), nu, to_var(L_S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<false>(to_var(L_Y), nu, L_S));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<false>(
      L_Y, to_var(nu), to_var(L_S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<false>(L_Y, to_var(nu), L_S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<false>(L_Y, nu, to_var(L_S)));

  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<true>(
      to_var(L_Y), to_var(nu), to_var(L_S)));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<true>(
      to_var(L_Y), to_var(nu), L_S));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<true>(
      to_var(L_Y), nu, to_var(L_S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<true>(to_var(L_Y), nu, L_S));
  test::check_varis_on_stack(stan::math::inv_wishart_cholesky_lpdf<true>(
      L_Y, to_var(nu), to_var(L_S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<true>(L_Y, to_var(nu), L_S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_cholesky_lpdf<true>(L_Y, nu, to_var(L_S)));

  stan::math::recover_memory();
}
