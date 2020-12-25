#include <limits>
#include <vector>
#include <gtest/gtest.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/tools/promotion.hpp>

#include <stan/math/rev.hpp>
#include <stan/math/prim/fun/sign.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/log1m.hpp>

template <typename T1, typename T2, typename T3, typename T4>
inline typename boost::math::tools::promote_args<T1, T2, T3, T4>::type
skew_de_cdf_test(
    const T1& y, const T2& mu, const T3& sigma, const T4& tau) {
  using std::log;
  using std::exp;

  if (y < mu) {
    return log(tau) - 2/sigma * (1-tau)*(mu -y);
  }
  else {
    return log1m((1- tau) * exp(-2/sigma * tau * (y-mu)));
  }
}


TEST(ProbDistributionsSkewedDoubleExponential, lpdf_computes_correct_gradients) {
  using stan::math::skew_double_exponential_lcdf;

        stan::math::var y = .1;
        stan::math::var mu = 1.2;
        stan::math::var sigma = 12.1;
        stan::math::var tau = 0.23;

        stan::math::var lp = skew_double_exponential_lcdf(y, mu, sigma, tau);
        std::vector<stan::math::var> theta;
        theta.push_back(y);
        theta.push_back(mu);
        theta.push_back(sigma);
        theta.push_back(tau);
        std::vector<double> g;
        lp.grad(theta, g);

        std::cout << " f = " << skew_double_exponential_lcdf(y, mu, sigma, tau) << std::endl;
        std::cout << " d.f / d.y = " << g[0] << std::endl;
        std::cout << " d.f / d.mu = " << g[1] << std::endl;
        std::cout << " d.f / d.sig = " << g[2] << std::endl;
        std::cout << " d.f / d.tau = " << g[3] << std::endl;

}
