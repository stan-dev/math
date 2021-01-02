#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, lub_mat_constrain) {
  using stan::scalar_type_t;
  using stan::math::lub_constrain;
  using stan::math::promote_scalar_t;
  auto tester = [](const auto& x) { return lub_constrain(x, -3.0, 3.0); };

  auto promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    x_scalar lb = -3.0;
    x_scalar ub = 3.0;
    return lub_constrain(x, lb, ub);
  };

  auto multi_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using lub_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar lb = -3.0;
    x_scalar ub = 3.0;
    lub_mat_t lb_mat = lub_mat_t::Constant(x.rows(), x.cols(), lb);
    lub_mat_t ub_mat = lub_mat_t::Constant(x.rows(), x.cols(), ub);
    return lub_constrain(x, lb_mat, ub_mat);
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(tester, A);
  stan::test::expect_ad_matvar(tester, A);
  stan::test::expect_ad(promote_tester, A);
  stan::test::expect_ad_matvar(promote_tester, A);
  stan::test::expect_ad_matvar(multi_promote_tester, A);
}

TEST(mathMixMatFun, lub_lp_mat_constrain) {
  using stan::scalar_type_t;
  using stan::math::lub_constrain;
  using stan::math::promote_scalar_t;
  auto tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    x_scalar lp = 0;
    lub_constrain(x, -3.0, 3.0, lp);
    return lp;
  };

  auto promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    x_scalar lb = -3.0;
    x_scalar ub = 3.0;
    x_scalar lp = 0;
    lub_constrain(x, lb, ub, lp);
    return lp;
  };

  auto multi_promote_tester = [](const auto& x) {
    using x_scalar = scalar_type_t<decltype(x)>;
    using lub_mat_t = Eigen::Matrix<x_scalar, -1, -1>;
    x_scalar lb = -3.0;
    x_scalar ub = 3.0;
    x_scalar lp = 0;
    lub_mat_t lb_mat = lub_mat_t::Constant(x.rows(), x.cols(), lb);
    lub_mat_t ub_mat = lub_mat_t::Constant(x.rows(), x.cols(), ub);
    lub_constrain(x, lb_mat, ub_mat, lp);
    return lp;
  };

  Eigen::MatrixXd A(Eigen::MatrixXd::Random(2, 2));
  stan::test::expect_ad(tester, A);
  stan::test::expect_ad_matvar(tester, A);
  stan::test::expect_ad(promote_tester, A);
  stan::test::expect_ad_matvar(promote_tester, A);
  stan::test::expect_ad_matvar(multi_promote_tester, A);
}
