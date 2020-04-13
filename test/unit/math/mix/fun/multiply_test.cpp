#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, multiply) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::multiply(x, y); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-1;
  tols.hessian_fvar_hessian_ = 1e-1;

  double a = 10;
  Eigen::VectorXd v1(1);
  v1 << 3;
  Eigen::RowVectorXd rv1(1);
  rv1 << -2;
  Eigen::MatrixXd m11(1, 1);
  m11 << 1.5;
  stan::test::expect_ad(tols, f, a, v1);
  stan::test::expect_ad(tols, f, v1, a);
  stan::test::expect_ad(tols, f, a, rv1);
  stan::test::expect_ad(tols, f, rv1, a);
  stan::test::expect_ad(tols, f, rv1, v1);
  stan::test::expect_ad(tols, f, v1, rv1);
  stan::test::expect_ad(tols, f, a, m11);
  stan::test::expect_ad(tols, f, m11, a);
  stan::test::expect_ad(tols, f, m11, v1);
  stan::test::expect_ad(tols, f, rv1, m11);
  stan::test::expect_ad(tols, f, m11, m11);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, a, v0);
  stan::test::expect_ad(f, v0, a);
  stan::test::expect_ad(f, a, rv0);
  stan::test::expect_ad(f, rv0, a);
  stan::test::expect_ad(f, a, m00);
  stan::test::expect_ad(f, m00, a);
  stan::test::expect_ad(f, m00, v0);
  stan::test::expect_ad(f, rv0, v0);
  stan::test::expect_ad(f, v0, rv0);
  stan::test::expect_ad(f, rv0, m00);
  stan::test::expect_ad(f, m00, m00);

  Eigen::VectorXd v(2);
  v << 100, -3;
  Eigen::RowVectorXd rv(2);
  rv << 100, -3;
  Eigen::MatrixXd m(2, 2);
  m << 100, 0, -3, 4;
  stan::test::expect_ad(tols, f, a, v);
  stan::test::expect_ad(tols, f, v, a);
  stan::test::expect_ad(tols, f, a, rv);
  stan::test::expect_ad(tols, f, rv, a);
  stan::test::expect_ad(tols, f, rv, v);
  stan::test::expect_ad(tols, f, v, rv);
  stan::test::expect_ad(tols, f, a, m);
  stan::test::expect_ad(tols, f, m, a);
  stan::test::expect_ad(tols, f, m, v);
  stan::test::expect_ad(tols, f, rv, m);
  stan::test::expect_ad(tols, f, m, m);

  Eigen::RowVectorXd d1(3);
  d1 << 1, 3, -5;
  Eigen::VectorXd d2(3);
  d2 << 4, -2, -1;
  stan::test::expect_ad(tols, f, d1, d2);
  stan::test::expect_ad(tols, f, d2, d1);

  Eigen::MatrixXd u(3, 2);
  u << 1, 3, -5, 4, -2, -1;
  Eigen::MatrixXd u_tr = u.transpose();
  Eigen::VectorXd vv(2);
  vv << -2, 4;
  Eigen::RowVectorXd rvv(3);
  rvv << -2, 4, 1;
  stan::test::expect_ad(tols, f, u, u_tr);
  stan::test::expect_ad(tols, f, u_tr, u);
  stan::test::expect_ad(tols, f, u, vv);
  stan::test::expect_ad(tols, f, rvv, u);

  // exception cases
  // can't compile mismatched dimensions, so no tests
}

template <typename T>
void instantiate_multiply() {
  Eigen::Matrix<double, -1, -1> d(2, 2);
  d << 1, 2, 3, 4;
  Eigen::Matrix<T, -1, -1> v(2, 2);
  v << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<double>, -1, -1> cd(2, 2);
  cd << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<T>, -1, -1> cv(2, 2);
  cv << 1, 2, 3, 4;

  auto d_d = d * d;
  auto d_v = d * v;
  auto d_cd = d * cd;
  auto d_cv = d * cv;

  auto v_d = v * d;
  auto v_v = v * v;
  auto v_cd = v * cd;
  auto v_cv = v * cv;

  auto cd_d = cd * d;
  auto cd_v = cd * v;
  auto cd_cd = cd * cd;
  auto cd_cv = cd * cv;

  auto cv_d = cv * d;
  auto cv_v = cv * v;
  auto cv_cd = cv * cd;
  auto cv_cv = cv * cv;
}

TEST(mathMix, multiplicationPatterns) {
  using stan::math::fvar;
  using stan::math::var;
  instantiate_multiply<double>();
  instantiate_multiply<var>();
  instantiate_multiply<fvar<double>>();
  instantiate_multiply<fvar<fvar<double>>>();
  instantiate_multiply<fvar<var>>();
  instantiate_multiply<fvar<fvar<var>>>();
}
