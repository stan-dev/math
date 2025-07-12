
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/tools/promotion.hpp>
#include <limits>
#include <vector>
#include <gtest/gtest.h>

namespace {
template <typename T1, typename T2, typename T3, typename T4>
inline auto skew_de_ccdf_test(const T1& y, const T2& mu, const T3& sigma,
                              const T4& tau) {
  using stan::math::log1m;
  using stan::math::log1m_exp;
  using std::exp;
  using std::log;

  if (y < mu) {
    return stan::math::log1m_exp(stan::math::log(tau)
                                 - 2.0 / sigma * (1.0 - tau) * (mu - y));
  } else {
    return stan::math::log1m_exp(stan::math::log1m(
        (1 - tau) * stan::math::exp(-2 / sigma * tau * (y - mu))));
  }
}

TEST(RevProbDistributionsSkewedDoubleExponential,
     lccdf_computes_correct_gradients) {
  using stan::math::skew_double_exponential_lccdf;

  for (double ys : {-1.7, 0.2, 0.5, 0.9, 1.1, 3.2, 8.3}) {
    for (double mus : {-1.8, 0.1, 0.55, 0.89, 1.3, 4.2, 9.3}) {
      for (double sigmas : {0.1, 0.5, 1.1, 10.1}) {
        for (double taus : {0.01, 0.1, 0.5, 0.6, 0.9, 0.99}) {
          stan::math::var y = ys;
          stan::math::var mu = mus;
          stan::math::var sigma = sigmas;
          stan::math::var tau = taus;

          stan::math::var lp = skew_double_exponential_lccdf(y, mu, sigma, tau);
          std::vector<stan::math::var> theta;
          theta.push_back(y);
          theta.push_back(mu);
          theta.push_back(sigma);
          theta.push_back(tau);
          std::vector<double> grads;
          lp.grad(theta, grads);

          stan::math::var y_true = ys;
          stan::math::var mu_true = mus;
          stan::math::var sigma_true = sigmas;
          stan::math::var tau_true = taus;

          stan::math::var lp_test
              = skew_de_ccdf_test(y_true, mu_true, sigma_true, tau_true);
          std::vector<stan::math::var> theta_true;
          theta_true.push_back(y_true);
          theta_true.push_back(mu_true);
          theta_true.push_back(sigma_true);
          theta_true.push_back(tau_true);
          std::vector<double> grads_true;
          lp_test.grad(theta_true, grads_true);

          EXPECT_NEAR(grads_true[0], grads[0], 0.01);
          EXPECT_NEAR(grads_true[1], grads[1], 0.01);
          EXPECT_NEAR(grads_true[2], grads[2], 0.01);
          EXPECT_NEAR(grads_true[3], grads[3], 0.01);
        }
      }
    }
  }
}
}  // namespace
