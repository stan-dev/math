#include <stan/math/mix.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradMixMath, ffv_gp_periodic_cov1_eigen_rvec_nan_values) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<Eigen::Matrix<double, 1, -1> > x_vec_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x_vec_fv(5);
  for (size_t i = 0; i < x_vec_d.size(); ++i) {
    x_vec_d[i].resize(1, 5);
    x_vec_fv[i].resize(1, 5);
    x_vec_d[i] << 1, 2, 3, 4, 5;
    x_vec_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x_vec_d_bad(x_vec_d);
  x_vec_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x_vec_fv_bad(x_vec_fv);
  x_vec_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv, sigma_d_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x_vec_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x_vec_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_vec_invalid_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > a(1.0, 1.0);
  fvar<fvar<var> > b(2.0, 1.0);
  fvar<fvar<var> > c(3.0, 1.0);
  fvar<fvar<var> > d(4.0, 1.0);
  fvar<fvar<var> > e(0.0, 1.0);

  std::vector<fvar<fvar<var> > > x1_fv;
  x1_fv.push_back(e);
  x1_fv.push_back(a);
  x1_fv.push_back(b);
  x1_fv.push_back(c);
  x1_fv.push_back(d);

  std::vector<fvar<fvar<var> > > x2_fv;
  x2_fv.push_back(2.);
  x2_fv.push_back(3.);
  x2_fv.push_back(1.);
  x2_fv.push_back(6.);
  x2_fv.push_back(-1.);

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x1_d;
  x1_d.push_back(1);
  x1_d.push_back(2);
  x1_d.push_back(3);
  x1_d.push_back(4);
  x1_d.push_back(5);

  std::vector<double> x2_d;
  x2_d.push_back(1);
  x2_d.push_back(2);
  x2_d.push_back(3);
  x2_d.push_back(4);
  x2_d.push_back(5);

  std::vector<double> x1_d_bad(x1_d);
  x1_d_bad[1] = nan_d;

  std::vector<double> x2_d_bad(x2_d);
  x2_d_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x1_fv_bad(x1_fv);
  x1_fv_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x2_fv_bad(x2_fv);
  x2_fv_bad[1] = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_vec_invalid_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > a(1.0, 1.0);
  fvar<fvar<var> > b(2.0, 1.0);
  fvar<fvar<var> > c(3.0, 1.0);
  fvar<fvar<var> > d(4.0, 1.0);
  fvar<fvar<var> > e(0.0, 1.0);

  std::vector<fvar<fvar<var> > > x1_fv;
  x1_fv.push_back(e);
  x1_fv.push_back(a);
  x1_fv.push_back(b);
  x1_fv.push_back(c);
  x1_fv.push_back(d);

  std::vector<fvar<fvar<var> > > x2_fv;
  x2_fv.push_back(2.);
  x2_fv.push_back(3.);
  x2_fv.push_back(1.);
  x2_fv.push_back(6.);
  x2_fv.push_back(-1.);

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<double> x1_d;
  x1_d.push_back(1);
  x1_d.push_back(2);
  x1_d.push_back(3);
  x1_d.push_back(4);
  x1_d.push_back(5);

  std::vector<double> x2_d;
  x2_d.push_back(1);
  x2_d.push_back(2);
  x2_d.push_back(3);
  x2_d.push_back(4);
  x2_d.push_back(5);

  std::vector<double> x1_d_bad(x1_d);
  x1_d_bad[1] = nan_d;

  std::vector<double> x2_d_bad(x2_d);
  x2_d_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x1_fv_bad(x1_fv);
  x1_fv_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x2_fv_bad(x2_fv);
  x2_fv_bad[1] = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_vec_nan_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > a(1.0, 1.0);
  fvar<fvar<var> > b(2.0, 1.0);
  fvar<fvar<var> > c(3.0, 1.0);
  fvar<fvar<var> > d(4.0, 1.0);
  fvar<fvar<var> > e(0.0, 1.0);

  std::vector<fvar<fvar<var> > > x1_fv;
  x1_fv.push_back(e);
  x1_fv.push_back(a);
  x1_fv.push_back(b);
  x1_fv.push_back(c);
  x1_fv.push_back(d);

  std::vector<fvar<fvar<var> > > x2_fv;
  x2_fv.push_back(2.);
  x2_fv.push_back(3.);
  x2_fv.push_back(1.);
  x2_fv.push_back(6.);
  x2_fv.push_back(-1.);

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<double> x1_d;
  x1_d.push_back(1);
  x1_d.push_back(2);
  x1_d.push_back(3);
  x1_d.push_back(4);
  x1_d.push_back(5);

  std::vector<double> x2_d;
  x2_d.push_back(1);
  x2_d.push_back(2);
  x2_d.push_back(3);
  x2_d.push_back(4);
  x2_d.push_back(5);

  std::vector<double> x1_d_bad(x1_d);
  x1_d_bad[1] = nan_d;

  std::vector<double> x2_d_bad(x2_d);
  x2_d_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x1_fv_bad(x1_fv);
  x1_fv_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x2_fv_bad(x2_fv);
  x2_fv_bad[1] = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_vec_nan_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > a(1.0, 1.0);
  fvar<fvar<var> > b(2.0, 1.0);
  fvar<fvar<var> > c(3.0, 1.0);
  fvar<fvar<var> > d(4.0, 1.0);
  fvar<fvar<var> > e(0.0, 1.0);

  std::vector<fvar<fvar<var> > > x1_fv;
  x1_fv.push_back(e);
  x1_fv.push_back(a);
  x1_fv.push_back(b);
  x1_fv.push_back(c);
  x1_fv.push_back(d);

  std::vector<fvar<fvar<var> > > x2_fv;
  x2_fv.push_back(2.);
  x2_fv.push_back(3.);
  x2_fv.push_back(1.);
  x2_fv.push_back(6.);
  x2_fv.push_back(-1.);

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<double> x1_d;
  x1_d.push_back(1);
  x1_d.push_back(2);
  x1_d.push_back(3);
  x1_d.push_back(4);
  x1_d.push_back(5);

  std::vector<double> x2_d;
  x2_d.push_back(1);
  x2_d.push_back(2);
  x2_d.push_back(3);
  x2_d.push_back(4);
  x2_d.push_back(5);

  std::vector<double> x1_d_bad(x1_d);
  x1_d_bad[1] = nan_d;

  std::vector<double> x2_d_bad(x2_d);
  x2_d_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x1_fv_bad(x1_fv);
  x1_fv_bad[1] = nan_d;

  std::vector<fvar<fvar<var> > > x2_fv_bad(x2_fv);
  x2_fv_bad[1] = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_vec_invalid_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_fv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(5, 1);
    x2_fv[i].resize(5, 1);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_vec_invalid_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_fv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(5, 1);
    x2_fv[i].resize(5, 1);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_vec_nan_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_fv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(5, 1);
    x2_fv[i].resize(5, 1);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_vec_nan_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_fv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(5, 1);
    x2_fv[i].resize(5, 1);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_rvec_invalid_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_fv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 5);
    x2_fv[i].resize(1, 5);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_rvec_invalid_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  fvar<fvar<var> > sigma_bad(-1.0, 1.0);
  fvar<fvar<var> > l_bad(-1.0, 1.0);
  fvar<fvar<var> > p_bad(-0.7, 1.0);
  double sigma_d_bad(-2.0);
  double l_d_bad(-1.0);
  double p_d_bad(-0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_fv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 5);
    x2_fv[i].resize(1, 5);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_rvec_nan_values_1) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_fv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 5);
    x2_fv[i].resize(1, 5);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_rvec_nan_values_2) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma(2.0, 1.0);
  fvar<fvar<var> > l(1.0, 1.0);
  fvar<fvar<var> > p(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  double nan_d = std::numeric_limits<double>::quiet_NaN();

  fvar<fvar<var> > sigma_bad(nan_d, 1.0);
  fvar<fvar<var> > l_bad(nan_d, 1.0);
  fvar<fvar<var> > p_bad(nan_d, 1.0);
  double sigma_d_bad(nan_d);
  double l_d_bad(nan_d);
  double p_d_bad(nan_d);

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_fv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 5);
    x2_fv[i].resize(1, 5);
    x2_d[i] << 5, 2, 1, 4, 5;
    x2_fv[i] << 4, 2, 2, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d_bad(x1_d);
  x1_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_fv_bad(x1_fv);
  x1_fv_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d_bad(x2_d);
  x2_d_bad[1](1) = nan_d;

  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_fv_bad(x2_fv);
  x2_fv_bad[1](1) = nan_d;

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_bad, l_d_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_fv_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_fv_bad, sigma_d_bad, l_d_bad,
                                     p_d_bad),
               std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d, l_d, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_fv_bad, x2_d_bad, sigma_d_bad, l_d_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_d_bad, l_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d, l, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_d_bad, l_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d, p_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_d_bad, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d, p_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_d_bad, p_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l, p_d_bad),
      std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_fv_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);

  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_d_bad, sigma_bad, l_bad, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d),
               std::domain_error);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l, p_d_bad),
               std::domain_error);
  EXPECT_THROW(
      cov = gp_periodic_cov(x1_d_bad, x2_d_bad, sigma_bad, l_bad, p_d_bad),
      std::domain_error);
}

TEST(AgradMixMath, fv_gp_periodic_cov2_eigen_vec_dim_error) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<var> sigma_fv(2.0, 1.0);
  fvar<var> l_fv(1.0, 1.0);
  fvar<var> p_fv(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<var>, -1, 1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_fv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<var>, -1, 1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(4, 1);
    x2_fv[i].resize(4, 1);
    x2_d[i] << 5, 2, 1, 4;
    x2_fv[i] << 4, 2, 2, 4;
  }

  Eigen::Matrix<fvar<var>, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_fv, p_fv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_fv, p_fv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_fv, p_fv),
               std::invalid_argument);
}

TEST(AgradMixMath, fv_gp_periodic_cov2_eigen_rvec_dim_error) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<var> sigma_fv(2.0, 1.0);
  fvar<var> l_fv(1.0, 1.0);
  fvar<var> p_fv(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<var>, 1, -1> > x1_fv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_fv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_fv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<var>, 1, -1> > x2_fv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 4);
    x2_fv[i].resize(1, 4);
    x2_d[i] << 5, 2, 1, 4;
    x2_fv[i] << 4, 2, 2, 4;
  }

  Eigen::Matrix<fvar<var>, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_fv, sigma_fv, l_fv, p_fv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_fv, sigma_fv, l_fv, p_fv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_fv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_d, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_d, l_fv, p_fv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_fv, x2_d, sigma_fv, l_fv, p_fv),
               std::invalid_argument);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_vec_dim_error) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma_ffv(2.0, 1.0);
  fvar<fvar<var> > l_ffv(1.0, 1.0);
  fvar<fvar<var> > p_ffv(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  std::vector<Eigen::Matrix<double, -1, 1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x1_ffv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(5, 1);
    x1_ffv[i].resize(5, 1);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_ffv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, -1, 1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, -1, 1> > x2_ffv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(4, 1);
    x2_ffv[i].resize(4, 1);
    x2_d[i] << 5, 2, 1, 4;
    x2_ffv[i] << 4, 2, 2, 4;
  }

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);
}

TEST(AgradMixMath, ffv_gp_periodic_cov2_eigen_rvec_dim_error) {
  using stan::math::fvar;
  using stan::math::gp_periodic_cov;
  using stan::math::var;

  fvar<fvar<var> > sigma_ffv(2.0, 1.0);
  fvar<fvar<var> > l_ffv(1.0, 1.0);
  fvar<fvar<var> > p_ffv(0.7, 1.0);
  double sigma_d(2.0);
  double l_d(1.0);
  double p_d(0.7);

  std::vector<Eigen::Matrix<double, 1, -1> > x1_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x1_ffv(5);
  for (size_t i = 0; i < x1_d.size(); ++i) {
    x1_d[i].resize(1, 5);
    x1_ffv[i].resize(1, 5);
    x1_d[i] << 1, 2, 3, 4, 5;
    x1_ffv[i] << 1, 2, 3, 4, 5;
  }

  std::vector<Eigen::Matrix<double, 1, -1> > x2_d(5);
  std::vector<Eigen::Matrix<fvar<fvar<var> >, 1, -1> > x2_ffv(5);
  for (size_t i = 0; i < x2_d.size(); ++i) {
    x2_d[i].resize(1, 4);
    x2_ffv[i].resize(1, 4);
    x2_d[i] << 5, 2, 1, 4;
    x2_ffv[i] << 4, 2, 2, 4;
  }

  Eigen::Matrix<fvar<fvar<var> >, -1, -1> cov;
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_d, x2_ffv, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_ffv, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);

  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_d, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_ffv, p_d),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_d, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_d, l_ffv, p_ffv),
               std::invalid_argument);
  EXPECT_THROW(cov = gp_periodic_cov(x1_ffv, x2_d, sigma_ffv, l_ffv, p_ffv),
               std::invalid_argument);
}
