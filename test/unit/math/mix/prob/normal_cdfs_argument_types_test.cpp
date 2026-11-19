#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <vector>

namespace normal_cdfs_types_test {
// values of the double overloads
inline double std_cdf(const Eigen::VectorXd& y) {
  return stan::math::std_normal_cdf(y);
}
inline double norm_cdf(const Eigen::VectorXd& y, double mu, double sigma) {
  return stan::math::normal_cdf(y, mu, sigma);
}
inline double lnorm_cdf(const Eigen::VectorXd& y, double mu, double sigma) {
  return stan::math::lognormal_cdf(y, mu, sigma);
}
}  // namespace normal_cdfs_types_test

TEST_F(AgradRev, mathMixScalFun_normal_cdfs_nested_fvar) {
  using stan::math::fvar;
  using stan::math::var;
  using ffv = fvar<fvar<var>>;
  ffv y = 0.3;
  y.d_.val_ = 1.0;
  ffv sigma = 1.2;
  ffv p = stan::math::std_normal_cdf(y);
  EXPECT_NEAR(stan::math::Phi(0.3), p.val_.val_.val(), 1e-15);
  EXPECT_NEAR(std::exp(stan::math::std_normal_lpdf(0.3)), p.d_.val_.val(),
              1e-15);
  ffv q = stan::math::normal_cdf(y, 0.0, sigma);
  EXPECT_NEAR(stan::math::Phi(0.3 / 1.2), q.val_.val_.val(), 1e-15);
  EXPECT_NEAR(std::exp(stan::math::normal_lpdf(0.3, 0.0, 1.2)), q.d_.val_.val(),
              1e-15);
  ffv r = stan::math::lognormal_cdf(sigma, y, sigma);
  EXPECT_NEAR(stan::math::lognormal_cdf(1.2, 0.3, 1.2), r.val_.val_.val(),
              1e-15);
}

TEST_F(AgradRev, mathMixScalFun_normal_cdfs_std_vector_var) {
  using stan::math::var;
  Eigen::VectorXd yd(2);
  yd << 0.1, 0.2;
  std::vector<var> y{0.1, 0.2};
  var mu = 0.0, sigma = 1.0;
  EXPECT_NEAR(normal_cdfs_types_test::std_cdf(yd),
              stan::math::std_normal_cdf(y).val(), 1e-15);
  EXPECT_NEAR(normal_cdfs_types_test::norm_cdf(yd, 0.0, 1.0),
              stan::math::normal_cdf(y, mu, sigma).val(), 1e-15);
  EXPECT_NEAR(normal_cdfs_types_test::lnorm_cdf(yd, 0.0, 1.0),
              stan::math::lognormal_cdf(y, mu, sigma).val(), 1e-15);
  // vector in the location or the scale only
  EXPECT_NEAR(stan::math::normal_cdf(0.5, std::vector<double>{0.1, 0.2}, 1.0),
              stan::math::normal_cdf(0.5, y, sigma).val(), 1e-15);
  EXPECT_NEAR(
      stan::math::lognormal_cdf(1.5, 0.0, std::vector<double>{0.1, 0.2}),
      stan::math::lognormal_cdf(1.5, mu, y).val(), 1e-15);
}

TEST_F(AgradRev, mathMixScalFun_normal_cdfs_eigen_fvar) {
  using stan::math::fvar;
  Eigen::VectorXd yd(2);
  yd << 0.1, 0.2;
  Eigen::Matrix<fvar<double>, -1, 1> y(2);
  y << 0.1, 0.2;
  EXPECT_NEAR(normal_cdfs_types_test::std_cdf(yd),
              stan::math::std_normal_cdf(y).val_, 1e-15);
  EXPECT_NEAR(normal_cdfs_types_test::norm_cdf(yd, 0.0, 1.0),
              stan::math::normal_cdf(y, 0.0, 1.0).val_, 1e-15);
  // the product is about 7e-147, so compare relative to it
  EXPECT_NEAR(1.0,
              stan::math::lognormal_cdf(y, y, y).val_
                  / stan::math::lognormal_cdf(yd, yd, yd),
              1e-14);
}

TEST_F(AgradRev, mathMixScalFun_normal_cdfs_varmat) {
  using stan::math::var;
  using stan::math::var_value;
  Eigen::VectorXd yd = Eigen::VectorXd::Constant(2, 0.5);
  yd(1) = 1.5;

  // each function: var_value<VectorXd> against Matrix<var>
  auto check = [&](auto f) {
    Eigen::Matrix<var, -1, 1> y_mat = yd;
    var mu_m = 0.2, sigma_m = 1.3;
    var lp_m = f(y_mat, mu_m, sigma_m);
    lp_m.grad();
    const Eigen::VectorXd y_adj = y_mat.adj();
    const double mu_adj = mu_m.adj(), sigma_adj = sigma_m.adj();
    stan::math::set_zero_all_adjoints();
    var_value<Eigen::VectorXd> y_var(yd);
    var mu = 0.2, sigma = 1.3;
    var lp = f(y_var, mu, sigma);
    lp.grad();
    EXPECT_NEAR(lp_m.val(), lp.val(), 1e-15);
    EXPECT_MATRIX_NEAR(y_adj, y_var.adj(), 1e-14);
    EXPECT_NEAR(mu_adj, mu.adj(), 1e-14);
    EXPECT_NEAR(sigma_adj, sigma.adj(), 1e-14);
    stan::math::set_zero_all_adjoints();
  };
  check([](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::normal_cdf(y, mu, sigma);
  });
  check([](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::lognormal_cdf(y, mu, sigma);
  });
  check([](const auto& y, const auto& mu, const auto& sigma) {
    return stan::math::std_normal_cdf(y) + 0 * (mu + sigma);
  });

  // row vectors and varmat in the location or the scale
  var_value<Eigen::RowVectorXd> vr(Eigen::RowVectorXd::Constant(2, 0.5));
  var_value<Eigen::VectorXd> vm(yd);
  var m = 0.0, s = 1.0;
  EXPECT_NEAR(
      stan::math::normal_cdf(0.0, 0.0, Eigen::RowVectorXd::Constant(2, 0.5)),
      stan::math::normal_cdf(0.0, 0.0, vr).val(), 1e-15);
  EXPECT_NEAR(
      stan::math::lognormal_cdf(Eigen::RowVectorXd::Constant(2, 0.5), 0.0,
                                Eigen::RowVectorXd::Constant(2, 0.5)),
      stan::math::lognormal_cdf(vr, m, vr).val(), 1e-15);
  EXPECT_NEAR(stan::math::lognormal_cdf(yd, yd, 1.0),
              stan::math::lognormal_cdf(vm, vm, s).val(), 1e-15);
}
