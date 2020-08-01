#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsSkewedDoubleExponential, test) {
  using stan::math::skew_double_exponential_lpdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  stan::math::var x(1);
  stan::math::var y(1);
  skew_double_exponential_lpdf(1, 0, x, y);
}
