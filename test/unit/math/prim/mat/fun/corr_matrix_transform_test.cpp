#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

using Eigen::Matrix;
using Eigen::Dynamic;

TEST(prob_transform,corr_matrix_j) {
  size_t K = 4;
  size_t K_choose_2 = 6; 
  Matrix<double,Dynamic,1> x(K_choose_2);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5;
  double lp = -12.9;
  Matrix<double,Dynamic,Dynamic> y = stan::math::corr_matrix_constrain(x,K,lp);
  Matrix<double,Dynamic,1> xrt = stan::math::corr_matrix_free(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}
TEST(prob_transform,corr_matrix_j2x2) {
  // tests K=2 boundary case, which has a different implementation
  size_t K = 2;
  size_t K_choose_2 = 1; 
  Matrix<double,Dynamic,1> x(K_choose_2);
  x << -1.3;
  double lp = -12.9;
  Matrix<double,Dynamic,Dynamic> y = stan::math::corr_matrix_constrain(x,K,lp);
  Matrix<double,Dynamic,1> xrt = stan::math::corr_matrix_free(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}

TEST(prob_transform,corr_matrix_constrain_exception) {
  unsigned int K = 4;
  unsigned int K_choose_2 = 6; 
  Matrix<double,Dynamic,1> x(K_choose_2-1);
  double lp = -12.9;

  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K), std::invalid_argument);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K, lp), std::invalid_argument);
  
  x.resize(K_choose_2+1);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K), std::invalid_argument);
  EXPECT_THROW(stan::math::corr_matrix_constrain(x, K, lp), std::invalid_argument);
}
TEST(prob_transform,corr_matrix_rt) {
  unsigned int K = 4;
  unsigned int K_choose_2 = 6; 
  Matrix<double,Dynamic,1> x(K_choose_2);
  x << -1.0, 2.0, 0.0, 1.0, 3.0, -1.5;
  Matrix<double,Dynamic,Dynamic> y = stan::math::corr_matrix_constrain(x,K);
  Matrix<double,Dynamic,1> xrt = stan::math::corr_matrix_free(y);
  EXPECT_EQ(x.size(), xrt.size());
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_FLOAT_EQ(x[i], xrt[i]);
  }
}
TEST(prob_transform,corr_matrix_free_exception) {
  Matrix<double,Dynamic,Dynamic> y;

  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(0,10);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(10,0);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);
  y.resize(1,2);
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::invalid_argument);

  y.resize(2,2);
  y << 0, 0, 0, 0;
  EXPECT_THROW(stan::math::corr_matrix_free(y), std::domain_error);
}
