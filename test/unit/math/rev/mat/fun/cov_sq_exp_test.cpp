#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <stan/math/prim/mat/fun/cov_sq_exp.hpp>
#include <vector>

TEST(MathPrimMat, cov_sq_exp_rev) {

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> grad;
  std::vector<stan::math::var> params;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {

      Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> x(3,1);
      x(0) = -2;
      x(1) = -1;
      x(2) = -0.5;
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x(i));
      params.push_back(x(j));      
 
      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));
      cov(i, j).grad(params, grad);

      double distance = x(i).val() - x(j).val();
      double sq_l = l.val() * l.val();
      double sq_distance = pow(distance, 2);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 
                        sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 
                        -distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val * 
                        distance/sq_l, grad[3])
        << "index: (" << i << ", " << j << ")";


      params.clear();
      grad.clear();
      stan::math::recover_memory();
    }
  }
}
