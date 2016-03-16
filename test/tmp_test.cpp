#include <stan/math/fwd/arr.hpp>
#include <gtest/gtest.h>

TEST(foo, 000) {
  int n = 5;
  double mu = 2;
  double phi = 1;
  EXPECT_FLOAT_EQ(-2.4327908, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}

TEST(foo, 100) {
  std::vector<int> n;
  n.push_back(5);
  double mu = 2;
  double phi = 1;
  EXPECT_FLOAT_EQ(-2.4327908, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}

TEST(foo, 010) {
  int n = 5;
  std::vector<double> mu;
  mu.push_back(2);
  double phi = 1;
  EXPECT_FLOAT_EQ(-2.4327908, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}


TEST(foo, 001) {
  int n = 5;
  double mu = 2;
  std::vector<double> phi;
  phi.push_back(1);
  EXPECT_FLOAT_EQ(-2.4327908, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}

TEST(foo, 001_zero_length) {
  int n = 5;
  double mu = 0;
  std::vector<double> phi;
  EXPECT_FLOAT_EQ(0, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}


TEST(foo, 111_zero_length) {
  std::vector<int> n;
  std::vector<double> mu;
  std::vector<double> phi;
  EXPECT_FLOAT_EQ(0, stan::math::neg_binomial_2_ccdf_log(n, mu, phi));
}



TEST(foo, fd_111_zero_length) {
  std::vector<int> n;
  std::vector<stan::math::fvar<double> > mu;
  std::vector<stan::math::fvar<double> > phi;
  EXPECT_FLOAT_EQ(1, stan::math::neg_binomial_2_cdf(n, mu, phi).val_);
}
