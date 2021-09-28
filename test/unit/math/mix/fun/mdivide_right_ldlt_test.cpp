#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, mdivideRightLdlt) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& x, const auto& y) {
    auto y_sym = stan::math::multiply(0.5, y + y.transpose()).eval();
    auto ldlt = stan::math::make_ldlt_factor(y_sym);
    return stan::math::mdivide_right_ldlt(x, ldlt);
  };

  Eigen::MatrixXd m00(0, 0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, rv0, m00);
  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, rv0, m00);

  Eigen::MatrixXd a(1, 1);
  a << 2;
  Eigen::MatrixXd b(1, 1);
  b << 3;
  stan::test::expect_ad(f, a, b);
  stan::test::expect_ad_matvar(f, a, b);

  Eigen::RowVectorXd g(1);
  g << 3;
  stan::test::expect_ad(f, g, a);
  stan::test::expect_ad_matvar(f, g, a);

  Eigen::MatrixXd c(2, 2);
  c << 2, 0.5, 0.5, 7;
  Eigen::MatrixXd d(2, 2);
  d << 2, 0, 0, 3;
  Eigen::MatrixXd ee(2, 2);
  ee << 12, 1, 1, 17;
  for (const auto& m1 : std::vector<Eigen::MatrixXd>{c, d, ee}) {
    for (const auto& m2 : std::vector<Eigen::MatrixXd>{c, d, ee}) {
      stan::test::expect_ad(f, m1, m2);
      stan::test::expect_ad_matvar(f, m1, m2);
    }
  }

  Eigen::VectorXd e(2);
  e << 2, 3;
  stan::test::expect_ad(f, e, c);
  stan::test::expect_ad_matvar(f, e, c);

  // ill-formed matrix inputs compile then throw at runtime
  Eigen::MatrixXd m33(3, 3);
  m33 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::MatrixXd m44(4, 4);
  m44 << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;

  Eigen::VectorXd v3(3);
  v3 << 1, 2, 3;

  Eigen::RowVectorXd rv3(3);
  rv3 << 1, 2, 3;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(2e-3, 2e-4);
  tols.hessian_fvar_hessian_ = relative_tolerance(2e-3, 2e-4);
  Eigen::RowVectorXd u(5);
  u << 62, 84, 84, 76, 108;
  Eigen::MatrixXd v(5, 5);
  v << 20, 8, -9, 7, 5, 8, 20, 0, 4, 4, -9, 0, 20, 2, 5, 7, 4, 2, 20, -5, 5, 4,
      5, -5, 20;
  stan::test::expect_ad(tols, f, u, v);
  stan::test::expect_ad_matvar(tols, f, u, v);

  // ill-formed inputs
  stan::test::expect_ad(f, m33, m44);         // wrong size
  stan::test::expect_ad(f, rv3, m44);         // wrong size
  stan::test::expect_ad(f, v3, m33);          // wrong type
  stan::test::expect_ad_matvar(f, m33, m44);  // wrong size
  stan::test::expect_ad_matvar(f, rv3, m44);  // wrong size
  stan::test::expect_ad_matvar(f, v3, m33);   // wrong type

  stan::math::recover_memory();
}
