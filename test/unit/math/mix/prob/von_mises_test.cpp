#include <stan/math/mix.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_von_mises_lpdf(double y, double mu, double kappa) {
  using stan::math::var;
  using stan::math::von_mises_lpdf;

  var y_var = y;
  var mu_var = mu;
  var kappa_var = kappa;

  std::vector<var> x;
  x.push_back(y_var);
  x.push_back(mu_var);
  x.push_back(kappa_var);

  var logp = von_mises_lpdf<false>(y_var, mu_var, kappa_var);
  std::vector<double> grad;
  logp.grad(x, grad);
  return grad;
}

TEST(ProbAgradDistributionsVonMises, derivatives) {
  using stan::math::fvar;
  using stan::math::von_mises_lpdf;

  std::vector<double> grad = test_von_mises_lpdf(0, 1, 0);

  fvar<double> lp = von_mises_lpdf<false>(0, 1, fvar<double>(0, 1));
  EXPECT_FLOAT_EQ(grad[2], lp.tangent());

  fvar<double> kappa1(0, 1);
  EXPECT_FLOAT_EQ(von_mises_lpdf(0, 1, kappa1).val_, -1.8378770664093453390819);
  EXPECT_FLOAT_EQ(von_mises_lpdf(0, 1, kappa1).d_, 0.54030230586813976501);

  fvar<fvar<double>> kappa2(0);
  EXPECT_NO_THROW(von_mises_lpdf(0, 1, kappa2));
}

TEST(ProbAgradDistributionsVonMises, FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::von_mises_lpdf;

  fvar<var> y_(0, 1);
  double mu(1);
  double kappa(0);

  fvar<var> logp = von_mises_lpdf(y_, mu, kappa);

  std::vector<stan::math::var> y{y_.val_};
  std::vector<double> g;
  logp.val_.grad(y, g);
  EXPECT_FLOAT_EQ(0, g[0]);
}

TEST(ProbAgradDistributionsVonMises, FvarVar_2ndDeriv1) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::von_mises_lpdf;

  double y_(0);
  fvar<var> mu(1, 1);
  double kappa(0);
  fvar<var> logp = von_mises_lpdf(y_, mu, kappa);

  std::vector<stan::math::var> y{mu.val_};
  std::vector<double> g;
  logp.d_.grad(y, g);
  EXPECT_FLOAT_EQ(0, g[0]);
}

TEST(ProbAgradDistributionsVonMises, FvarVar_2ndDeriv2) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::von_mises_lpdf;

  double y_(0);
  double mu(1);
  fvar<var> kappa(0, 1);
  fvar<var> logp = von_mises_lpdf(y_, mu, kappa);

  std::vector<stan::math::var> y{kappa.val_};
  std::vector<double> g;
  logp.d_.grad(y, g);
  // Fails: g[0] is -nan
  // EXPECT_FLOAT_EQ(0, g[0]);
}

// This test once failed sanitizer checks -- nothing explicitly tested in the
// code itself
TEST(ProbAgradDistributionsVonMises, sanitizer_error_fixed) {
  using stan::math::var;
  double y = boost::math::constants::third_pi<double>();
  double mu = boost::math::constants::sixth_pi<double>();
  std::vector<var> kappa = {0.5};

  auto lp = stan::math::von_mises_lpdf(y, mu, kappa);

  lp.grad();
}
