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

  stan::test::expect_ad_matvar(tols, f, a, v1);
  stan::test::expect_ad_matvar(tols, f, v1, a);
  stan::test::expect_ad_matvar(tols, f, a, rv1);
  stan::test::expect_ad_matvar(tols, f, rv1, a);
  stan::test::expect_ad_matvar(tols, f, rv1, v1);
  stan::test::expect_ad_matvar(tols, f, v1, rv1);
  stan::test::expect_ad_matvar(tols, f, a, m11);
  stan::test::expect_ad_matvar(tols, f, m11, a);
  stan::test::expect_ad_matvar(tols, f, m11, v1);
  stan::test::expect_ad_matvar(tols, f, rv1, m11);
  stan::test::expect_ad_matvar(tols, f, m11, m11);

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

  stan::test::expect_ad_matvar(f, a, v0);
  stan::test::expect_ad_matvar(f, v0, a);
  stan::test::expect_ad_matvar(f, a, rv0);
  stan::test::expect_ad_matvar(f, rv0, a);
  stan::test::expect_ad_matvar(f, a, m00);
  stan::test::expect_ad_matvar(f, m00, a);
  stan::test::expect_ad_matvar(f, m00, v0);
  stan::test::expect_ad_matvar(f, rv0, v0);
  stan::test::expect_ad_matvar(f, v0, rv0);
  stan::test::expect_ad_matvar(f, rv0, m00);
  stan::test::expect_ad_matvar(f, m00, m00);

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

  stan::test::expect_ad_matvar(tols, f, a, v);
  stan::test::expect_ad_matvar(tols, f, v, a);
  stan::test::expect_ad_matvar(tols, f, a, rv);
  stan::test::expect_ad_matvar(tols, f, rv, a);
  stan::test::expect_ad_matvar(tols, f, rv, v);
  stan::test::expect_ad_matvar(tols, f, v, rv);
  stan::test::expect_ad_matvar(tols, f, a, m);
  stan::test::expect_ad_matvar(tols, f, m, a);
  stan::test::expect_ad_matvar(tols, f, m, v);
  stan::test::expect_ad_matvar(tols, f, rv, m);
  stan::test::expect_ad_matvar(tols, f, m, m);

  Eigen::RowVectorXd d1(3);
  d1 << 1, 3, -5;
  Eigen::VectorXd d2(3);
  d2 << 4, -2, -1;
  stan::test::expect_ad(tols, f, d1, d2);
  stan::test::expect_ad(tols, f, d2, d1);

  stan::test::expect_ad_matvar(tols, f, d1, d2);
  stan::test::expect_ad_matvar(tols, f, d2, d1);

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

  stan::test::expect_ad_matvar(tols, f, u, u_tr);
  stan::test::expect_ad_matvar(tols, f, u_tr, u);
  stan::test::expect_ad_matvar(tols, f, u, vv);
  stan::test::expect_ad_matvar(tols, f, rvv, u);

  // exception cases
  // can't compile mismatched dimensions, so no tests
}

template <typename T>
void instantiate_multiply() {
  using stan::math::multiply;
  Eigen::Matrix<double, -1, -1> d_mat(2, 2);
  d_mat << 1, 2, 3, 4;
  Eigen::Matrix<T, -1, -1> v_mat(2, 2);
  v_mat << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<double>, -1, -1> cd_mat(2, 2);
  cd_mat << 1, 2, 3, 4;
  Eigen::Matrix<std::complex<T>, -1, -1> cv_mat(2, 2);
  cv_mat << 1, 2, 3, 4;

  auto d_d_mat = stan::math::eval(multiply(d_mat, d_mat));
  auto d_v_mat = stan::math::eval(multiply(d_mat, v_mat));
  auto d_cd_mat = stan::math::eval(multiply(d_mat, cd_mat));
  auto d_cv_mat = stan::math::eval(multiply(d_mat, cv_mat));

  auto v_d_mat = stan::math::eval(multiply(v_mat, d_mat));
  auto v_v_mat = stan::math::eval(multiply(v_mat, v_mat));
  auto v_cd_mat = stan::math::eval(multiply(v_mat, cd_mat));
  auto v_cv_mat = stan::math::eval(multiply(v_mat, cv_mat));

  auto cd_d_mat = stan::math::eval(multiply(cd_mat, d_mat));
  auto cd_v_mat = stan::math::eval(multiply(cd_mat, v_mat));
  auto cd_cd_mat = stan::math::eval(multiply(cd_mat, cd_mat));
  auto cd_cv_mat = stan::math::eval(multiply(cd_mat, cv_mat));

  auto cv_d_mat = stan::math::eval(multiply(cv_mat, d_mat));
  auto cv_v_mat = stan::math::eval(multiply(cv_mat, v_mat));
  auto cv_cd_mat = stan::math::eval(multiply(cv_mat, cd_mat));
  auto cv_cv_mat = stan::math::eval(multiply(cv_mat, cv_mat));

  Eigen::Matrix<double, -1, 1> d_vec(2);
  d_vec << 1, 2;
  Eigen::Matrix<T, -1, 1> v_vec(2);
  v_vec << 1, 2;
  Eigen::Matrix<std::complex<double>, -1, 1> cd_vec(2);
  cd_vec << 1, 2;
  Eigen::Matrix<std::complex<T>, -1, 1> cv_vec(2);
  cv_vec << 1, 2;

  auto d_d_vec_mat = stan::math::eval(multiply(d_mat, d_vec));
  auto d_v_vec_mat = stan::math::eval(multiply(d_mat, v_vec));
  auto d_cd_vec_mat = stan::math::eval(multiply(d_mat, cd_vec));
  auto d_cv_vec_mat = stan::math::eval(multiply(d_mat, cv_vec));

  auto v_d_vec_mat = stan::math::eval(multiply(v_mat, d_vec));
  auto v_v_vec_mat = stan::math::eval(multiply(v_mat, v_vec));
  auto v_cd_vec_mat = stan::math::eval(multiply(v_mat, cd_vec));
  auto v_cv_vec_mat = stan::math::eval(multiply(v_mat, cv_vec));

  auto cd_d_vec_mat = stan::math::eval(multiply(cd_mat, d_vec));
  auto cd_v_vec_mat = stan::math::eval(multiply(cd_mat, v_vec));
  auto cd_cd_vec_mat = stan::math::eval(multiply(cd_mat, cd_vec));
  auto cd_cv_vec_mat = stan::math::eval(multiply(cd_mat, cv_vec));

  auto cv_d_vec_mat = stan::math::eval(multiply(cv_mat, d_vec));
  auto cv_v_vec_mat = stan::math::eval(multiply(cv_mat, v_vec));
  auto cv_cd_vec_mat = stan::math::eval(multiply(cv_mat, cd_vec));
  auto cv_cv_vec_mat = stan::math::eval(multiply(cv_mat, cv_vec));

  Eigen::Matrix<double, 1, -1> d_rowvec(2);
  d_rowvec << 1, 2;
  Eigen::Matrix<T, 1, -1> v_rowvec(2);
  v_rowvec << 1, 2;
  Eigen::Matrix<std::complex<double>, 1, -1> cd_rowvec(2);
  cd_rowvec << 1, 2;
  Eigen::Matrix<std::complex<T>, 1, -1> cv_rowvec(2);
  cv_rowvec << 1, 2;

  auto d_d_dot_prod = stan::math::eval(multiply(d_vec, d_rowvec));
  auto d_v_dot_prod = stan::math::eval(multiply(d_vec, v_rowvec));
  auto d_cd_dot_prod = stan::math::eval(multiply(d_vec, cd_rowvec));
  auto d_cv_dot_prod = stan::math::eval(multiply(d_vec, cv_rowvec));

  auto v_d_dot_prod = stan::math::eval(multiply(v_vec, d_rowvec));
  auto v_v_dot_prod = stan::math::eval(multiply(v_vec, v_rowvec));
  auto v_cd_dot_prod = stan::math::eval(multiply(v_vec, cd_rowvec));
  auto v_cv_dot_prod = stan::math::eval(multiply(v_vec, cv_rowvec));

  auto cd_d_dot_prod = stan::math::eval(multiply(cd_vec, d_rowvec));
  auto cd_v_dot_prod = stan::math::eval(multiply(cd_vec, v_rowvec));
  auto cd_cd_dot_prod = stan::math::eval(multiply(cd_vec, cd_rowvec));
  auto cd_cv_dot_prod = stan::math::eval(multiply(cd_vec, cv_rowvec));

  auto cv_d_dot_prod = stan::math::eval(multiply(cv_vec, d_rowvec));
  auto cv_v_dot_prod = stan::math::eval(multiply(cv_vec, v_rowvec));
  auto cv_cd_dot_prod = stan::math::eval(multiply(cv_vec, cd_rowvec));
  auto cv_cv_dot_prod = stan::math::eval(multiply(cv_vec, cv_rowvec));

  auto d_d_outer_prod = stan::math::eval(multiply(d_rowvec, d_vec));
  auto d_v_outer_prod = stan::math::eval(multiply(d_rowvec, v_vec));
  auto d_cd_outer_prod = stan::math::eval(multiply(d_rowvec, cd_vec));
  auto d_cv_outer_prod = stan::math::eval(multiply(d_rowvec, cv_vec));

  auto v_d_outer_prod = stan::math::eval(multiply(v_rowvec, d_vec));
  auto v_v_outer_prod = stan::math::eval(multiply(v_rowvec, v_vec));
  auto v_cd_outer_prod = stan::math::eval(multiply(v_rowvec, cd_vec));
  auto v_cv_outer_prod = stan::math::eval(multiply(v_rowvec, cv_vec));

  auto cd_d_outer_prod = stan::math::eval(multiply(cd_rowvec, d_vec));
  auto cd_v_outer_prod = stan::math::eval(multiply(cd_rowvec, v_vec));
  auto cd_cd_outer_prod = stan::math::eval(multiply(cd_rowvec, cd_vec));
  auto cd_cv_outer_prod = stan::math::eval(multiply(cd_rowvec, cv_vec));

  auto cv_d_outer_prod = stan::math::eval(multiply(cv_rowvec, d_vec));
  auto cv_v_outer_prod = stan::math::eval(multiply(cv_rowvec, v_vec));
  auto cv_cd_outer_prod = stan::math::eval(multiply(cv_rowvec, cd_vec));
  auto cv_cv_outer_prod = stan::math::eval(multiply(cv_rowvec, cv_vec));

  auto d_d_rowvec_mat = stan::math::eval(multiply(d_rowvec, d_mat));
  auto d_v_rowvec_mat = stan::math::eval(multiply(d_rowvec, v_mat));
  auto d_cd_rowvec_mat = stan::math::eval(multiply(d_rowvec, cd_mat));
  auto d_cv_rowvec_mat = stan::math::eval(multiply(d_rowvec, cv_mat));

  auto v_d_rowvec_mat = stan::math::eval(multiply(v_rowvec, d_mat));
  auto v_v_rowvec_mat = stan::math::eval(multiply(v_rowvec, v_mat));
  auto v_cd_rowvec_mat = stan::math::eval(multiply(v_rowvec, cd_mat));
  auto v_cv_rowvec_mat = stan::math::eval(multiply(v_rowvec, cv_mat));

  auto cd_d_rowvec_mat = stan::math::eval(multiply(cd_rowvec, d_mat));
  auto cd_v_rowvec_mat = stan::math::eval(multiply(cd_rowvec, v_mat));
  auto cd_cd_rowvec_mat = stan::math::eval(multiply(cd_rowvec, cd_mat));
  auto cd_cv_rowvec_mat = stan::math::eval(multiply(cd_rowvec, cv_mat));

  auto cv_d_rowvec_mat = stan::math::eval(multiply(cv_rowvec, d_mat));
  auto cv_v_rowvec_mat = stan::math::eval(multiply(cv_rowvec, v_mat));
  auto cv_cd_rowvec_mat = stan::math::eval(multiply(cv_rowvec, cd_mat));
  auto cv_cv_rowvec_mat = stan::math::eval(multiply(cv_rowvec, cv_mat));

  auto d_d_mat_vec = stan::math::eval(multiply(d_mat, d_vec));
  auto d_v_mat_vec = stan::math::eval(multiply(d_mat, v_vec));
  auto d_cd_mat_vec = stan::math::eval(multiply(d_mat, cd_vec));
  auto d_cv_mat_vec = stan::math::eval(multiply(d_mat, cv_vec));

  auto v_d_mat_vec = stan::math::eval(multiply(v_mat, d_vec));
  auto v_v_mat_vec = stan::math::eval(multiply(v_mat, v_vec));
  auto v_cd_mat_vec = stan::math::eval(multiply(v_mat, cd_vec));
  auto v_cv_mat_vec = stan::math::eval(multiply(v_mat, cv_vec));

  auto cd_d_mat_vec = stan::math::eval(multiply(cd_mat, d_vec));
  auto cd_v_mat_vec = stan::math::eval(multiply(cd_mat, v_vec));
  auto cd_cd_mat_vec = stan::math::eval(multiply(cd_mat, cd_vec));
  auto cd_cv_mat_vec = stan::math::eval(multiply(cd_mat, cv_vec));

  auto cv_d_mat_vec = stan::math::eval(multiply(cv_mat, d_vec));
  auto cv_v_mat_vec = stan::math::eval(multiply(cv_mat, v_vec));
  auto cv_cd_mat_vec = stan::math::eval(multiply(cv_mat, cd_vec));
  auto cv_cv_mat_vec = stan::math::eval(multiply(cv_mat, cv_vec));
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
