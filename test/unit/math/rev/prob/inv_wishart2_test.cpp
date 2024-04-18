#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <gtest/gtest.h>
#include <string>

template <typename T_y, typename T_dof, typename T_scale>
void expect_propto_inv_wishart_lpdf(T_y W1, T_dof nu1, T_scale S1, T_y W2,
                                    T_dof nu2, T_scale S2,
                                    std::string message) {
  expect_eq_diffs(stan::math::inv_wishart_lpdf<false>(W1, nu1, S1),
                  stan::math::inv_wishart_lpdf<false>(W2, nu2, S2),
                  stan::math::inv_wishart_lpdf<true>(W1, nu1, S1),
                  stan::math::inv_wishart_lpdf<true>(W2, nu2, S2), message);
}

class AgradDistributionsInvWishart : public ::testing::Test {
 protected:
  virtual void SetUp() {
    Y1.resize(2, 2);
    Y1 << 2.011108, -11.20661, -11.20661, 112.94139;
    Y2.resize(2, 2);
    Y2 << 13.4, 12.2, 12.2, 11.5;

    nu1 = 3;
    nu2 = 5.3;

    S1.resize(2, 2);
    S1 << 1.848220, 1.899623, 1.899623, 12.751941;
    S2.resize(2, 2);
    S2 << 3.0, 1.4, 1.4, 7.0;
  }

  virtual void TearDown() { stan::math::recover_memory(); }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Y1;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Y2;
  double nu1;
  double nu2;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S1;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S2;
};

TEST_F(AgradDistributionsInvWishart, Propto) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(to_var(Y1), to_var(nu1), to_var(S1),
                                 to_var(Y2), to_var(nu2), to_var(S2),
                                 "var: y, nu, and sigma");
}
TEST_F(AgradDistributionsInvWishart, ProptoY) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(to_var(Y1), nu1, S1, to_var(Y2), nu1, S1,
                                 "var: y");
}
TEST_F(AgradDistributionsInvWishart, ProptoYNu) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(to_var(Y1), to_var(nu1), S1, to_var(Y2),
                                 to_var(nu2), S1, "var: y, and nu");
}
TEST_F(AgradDistributionsInvWishart, ProptoYSigma) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(to_var(Y1), nu1, to_var(S1), to_var(Y2), nu1,
                                 to_var(S2), "var: y and sigma");
}
TEST_F(AgradDistributionsInvWishart, ProptoNu) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(Y1, to_var(nu1), S1, Y1, to_var(nu2), S1,
                                 "var: nu");
}
TEST_F(AgradDistributionsInvWishart, ProptoNuSigma) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(Y1, to_var(nu1), to_var(S1), Y1, to_var(nu2),
                                 to_var(S2), "var: nu and sigma");
}
TEST_F(AgradDistributionsInvWishart, ProptoSigma) {
  using stan::math::to_var;
  expect_propto_inv_wishart_lpdf(Y1, nu1, to_var(S1), Y1, nu1, to_var(S2),
                                 "var: sigma");
}

TEST(InvWishart, check_varis_on_stack) {
  using stan::math::to_var;
  Eigen::MatrixXd W(2, 2);
  W << 2.011108, -11.20661, -11.20661, 112.94139;

  double nu = 3;

  Eigen::MatrixXd S(2, 2);
  S << 1.848220, 1.899623, 1.899623, 12.751941;

  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(to_var(W), to_var(nu), to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(to_var(W), to_var(nu), S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(to_var(W), nu, to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(to_var(W), nu, S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(W, to_var(nu), to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(W, to_var(nu), S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<false>(W, nu, to_var(S)));

  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(to_var(W), to_var(nu), to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(to_var(W), to_var(nu), S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(to_var(W), nu, to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(to_var(W), nu, S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(W, to_var(nu), to_var(S)));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(W, to_var(nu), S));
  test::check_varis_on_stack(
      stan::math::inv_wishart_lpdf<true>(W, nu, to_var(S)));

  stan::math::recover_memory();
}
