#include <stan/math/mix/scal.hpp>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_fun(double y) {
  using stan::math::std_normal_lpdf;
  using stan::math::var;

  var y_var = y;
  std::vector<var> x;
  x.push_back(y_var);

  var logp = std_normal_lpdf<false>(y_var);
  std::vector<double> grad;
  logp.grad(x, grad);
  return grad;
}

TEST(ProbAgradDistributionsStdNormal, derivatives) {
  using stan::math::fvar;
  using stan::math::std_normal_lpdf;

  std::vector<double> grad = test_fun(0);

  fvar<double> lp = std_normal_lpdf<false>(0);
  EXPECT_FLOAT_EQ(grad[2], lp.tangent());

  fvar<fvar<double>> y(1.0);
  fvar<double> x(1.0, 2.0);
  EXPECT_NO_THROW(std_normal_lpdf(y));
  EXPECT_FLOAT_EQ(std_normal_lpdf(x).val_, -1.418938533204672669541);
  EXPECT_FLOAT_EQ(std_normal_lpdf(x).d_, -2);
}

TEST(ProbAgradDistributionsStdNormal, FvarVar_1stDeriv) {
  using stan::math::fvar;
  using stan::math::std_normal_lpdf;
  using stan::math::var;

  fvar<var> y_(2, 1);
  fvar<var> logp = std_normal_lpdf(y_);

  AVEC y = createAVEC(y_.val_);
  VEC g;
  logp.val_.grad(y, g);
  EXPECT_FLOAT_EQ(-2, g[0]);
}
