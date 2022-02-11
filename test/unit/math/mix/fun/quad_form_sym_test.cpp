#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/ad_tolerances.hpp>

TEST(MathMixMatFun, quadFormSym) {
  auto f = [](const auto& x, const auto& y) {
    // symmetrize the input matrix
    auto x_sym = ((x + x.transpose()) * 0.5).eval();
    return stan::math::quad_form_sym(x_sym, y);
  };

  Eigen::MatrixXd a00;
  Eigen::MatrixXd a02(0, 2);
  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << -2;
  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  Eigen::MatrixXd b22(2, 2);
  b22 << -3, -2, -10, 112;
  Eigen::MatrixXd b23(2, 3);
  b23 << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b42(4, 2);
  b42 << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::MatrixXd a44(4, 4);
  a44 << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;

  Eigen::VectorXd v0(0);
  Eigen::VectorXd v1(1);
  v1 << 42;
  Eigen::VectorXd v2(2);
  v2 << -3, 13;
  Eigen::VectorXd v4(4);
  v4 << 100, 0, -3, 5;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 2e-1;
  tols.hessian_fvar_hessian_ = 2e-1;

  stan::test::expect_ad(f, a00, a00);
  stan::test::expect_ad(f, a00, a02);

  stan::test::expect_ad(f, a11, b11);
  stan::test::expect_ad(tols, f, a22, b22);
  stan::test::expect_ad(f, a22, b23);
  stan::test::expect_ad(tols, f, a44, b42);

  stan::test::expect_ad(f, a00, v0);
  stan::test::expect_ad(f, a11, v1);
  stan::test::expect_ad(f, a22, v2);
  stan::test::expect_ad(tols, f, a44, v4);

  // asymmetric case should throw

  auto g = [](const auto& x, const auto& y) {
    return stan::math::quad_form_sym(x, y);
  };

  stan::test::expect_ad(g, a02, a22);

  Eigen::MatrixXd u(4, 4);
  u << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;
  Eigen::MatrixXd v(4, 2);
  v << 100, 10, 0, 1, -3, -3, 5, 2;
  stan::test::expect_ad(g, u, v);
}

TEST(MathMixMatFun, quad_form_sym_2095) {
  Eigen::Matrix<stan::math::var, -1, -1> av(2, 2);
  Eigen::Matrix<stan::math::var, -1, -1> bv(2, 2);

  av << 1.25882993696386514, -0.03325949401909023, -0.03325949605426004,
      1.85523447220884385;

  bv << 2.55474619740069508, 0.66362927717720988, -1.92917223922387349,
      2.01853255256721731;

  Eigen::Matrix<stan::math::var, -1, -1> cv = stan::math::quad_form_sym(av, bv);
  EXPECT_FLOAT_EQ(0, cv(1, 0).val() - cv(0, 1).val());

  stan::math::recover_memory();

  //--------------------

  Eigen::MatrixXd ad(2, 2);
  Eigen::MatrixXd bd(2, 2);

  ad << 1.25882993696386514, -0.03325949401909023, -0.03325949605426004,
      1.85523447220884385;

  bd << 2.55474619740069508, 0.66362927717720988, -1.92917223922387349,
      2.01853255256721731;

  stan::math::quad_form_sym(ad, bd);

  auto f = [](const auto& x, const auto& y) {
    // symmetrize the input matrix;
    // expect_ad will perturb elements and cause it not to be symmetric
    auto x_sym = ((x + x.transpose()) * 0.5).eval();
    return stan::math::quad_form_sym(x_sym, y);
  };
  stan::test::expect_ad(f, ad, bd);
}
