#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <boost/math/tools/promotion.hpp>

template <typename T1, typename T2, typename T3, typename T4>
inline typename boost::math::tools::promote_args<T1, T2, T3, T4>::type blub(
    const T1& y, const T2& mu, const T3& sigma, const T4& tau) {
  using std::log;
  using std::pow;

  return log(2) + log(tau) + log(1 - tau) - log(sigma)
             - 2 * ((y < mu) ? (1 - tau) * (mu - y) : tau * (y - mu)) / sigma;
}


TEST(ProbDistributionsSkewedDoubleExponential, test) {
  using stan::math::skew_double_exponential_lpdf;

  stan::math::var y = 40.2;
  stan::math::var mu = 2.0;
  stan::math::var sigma = 2.0;
  stan::math::var tau = 0.23;

  stan::math::var lp = blub(y, mu, sigma, tau);
  std::vector<stan::math::var> theta;
  theta.push_back(y);
  theta.push_back(mu);
  theta.push_back(sigma);
  theta.push_back(tau);
  std::vector<double> g;
  lp.grad(theta, g);
  std::cout << " f = " << blub(y, mu, sigma, tau) << std::endl;
  std::cout << " d.f / d.y = " << g[0] << std::endl;
  std::cout << " d.f / d.mu = " << g[1] << std::endl;
  std::cout << " d.f / d.sig = " << g[2] << std::endl;
  std::cout << " d.f / d.tau = " << g[3] << std::endl;
}


TEST(ProbDistributionsSkewedDoubleExponential, test2) {
using stan::math::skew_double_exponential_lpdf;

  stan::math::var y2 = 40.2;
  stan::math::var mu2 = 2.0;
  stan::math::var sigma2 = 2.0;
  stan::math::var tau2 = 0.23;

  stan::math::var lp2 = skew_double_exponential_lpdf(y2, mu2, sigma2, tau2);
  std::vector<stan::math::var> theta2;
  theta2.push_back(y2);
  theta2.push_back(mu2);
  theta2.push_back(sigma2);
  theta2.push_back(tau2);
  std::vector<double> g2;
  lp2.grad(theta2, g2);
  std::cout << " f = " << skew_double_exponential_lpdf(
      y2, mu2, sigma2, tau2) << std::endl;
  std::cout << " d.f / d.y = " << g2[0] << std::endl;
  std::cout << " d.f / d.mu = " << g2[1] << std::endl;
  std::cout << " d.f / d.sig = " << g2[2] << std::endl;
  std::cout << " d.f / d.tau = " << g2[3] << std::endl;
}

