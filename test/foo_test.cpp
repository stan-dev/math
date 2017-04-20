#include <stan/math/mix/scal.hpp>
#include <gtest/gtest.h>

TEST(foo, bar) {
  double y = 0;
  double mu = 0;
  stan::math::fvar<stan::math::var> sigma = 1;


  stan::math::fvar<stan::math::var> lp = normal_lpdf(y, mu, sigma);
}
