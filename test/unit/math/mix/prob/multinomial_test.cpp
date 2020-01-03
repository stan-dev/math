#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

using Eigen::Dynamic;
using Eigen::Matrix;

TEST(ProbDistributionsMultinomial, fvar_var) {
  using stan::math::fvar;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<fvar<var>, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_log(ns, theta).val_.val());
  EXPECT_FLOAT_EQ(17.666666, stan::math::multinomial_log(ns, theta).d_.val());
}

TEST(ProbDistributionsMultinomial, fvar_fvar_var) {
  using stan::math::fvar;
  using stan::math::var;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<fvar<fvar<var> >, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  EXPECT_FLOAT_EQ(-2.002481,
                  stan::math::multinomial_log(ns, theta).val_.val_.val());
  EXPECT_FLOAT_EQ(17.666666,
                  stan::math::multinomial_log(ns, theta).d_.val_.val());
}
