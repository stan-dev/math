#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <vector>

TEST(ProbDistributionsMultinomial, fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<fvar<double>, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  for (int i = 0; i < 3; i++)
    theta(i).d_ = 1.0;

  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_lpmf(ns, theta).val_);
  EXPECT_FLOAT_EQ(17.666666, stan::math::multinomial_lpmf(ns, theta).d_);
}

TEST(ProbDistributionsMultinomial, fvar_fvar_double) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  std::vector<int> ns;
  ns.push_back(1);
  ns.push_back(2);
  ns.push_back(3);
  Matrix<fvar<fvar<double> >, Dynamic, 1> theta(3, 1);
  theta << 0.2, 0.3, 0.5;
  for (int i = 0; i < 3; i++)
    theta(i).d_.val_ = 1.0;

  EXPECT_FLOAT_EQ(-2.002481, stan::math::multinomial_lpmf(ns, theta).val_.val_);
  EXPECT_FLOAT_EQ(17.666666, stan::math::multinomial_lpmf(ns, theta).d_.val_);
}
