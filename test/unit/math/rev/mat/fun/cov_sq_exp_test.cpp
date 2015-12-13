#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <stan/math/prim/mat/fun/cov_sq_exp.hpp>
#include <vector>

TEST(RevMath, cov_sq_exp_vvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);      

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * -distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * distance/sq_l, grad[3])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var sigma = 0.2;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(x[i]);
      params.push_back(x[j]);      

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * -distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}


TEST(RevMath, cov_sq_exp_vdv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;

  double sigma = 0.2;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      stan::math::var l = 5;
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
      params.push_back(x[i]);
      params.push_back(x[j]);      

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * -distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * distance/sq_l, grad[2])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_vdd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  double sigma = 0.2;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      std::vector<stan::math::var> x(3);
      x[0] = -2;
      x[1] = -1;
      x[2] = -0.5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(x[i]);
      params.push_back(x[j]);      

      cov(i, j).grad(params, grad);

      double distance = x[i].val() - x[j].val();
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * -distance/sq_l, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * distance/sq_l, grad[1])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_dvv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
      params.push_back(l);
  
      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma.val() * sigma.val() * exp_val
                      * sq_distance/(sq_l * l.val()), grad[1])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_dvd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double l = 5;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var sigma = 0.2;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(sigma);
  
      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l);
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma.val()) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(2 * sigma.val() * exp_val, grad[0])
        << "index: (" << i << ", " << j << ")";

      stan::math::recover_memory();
    }
  }
}

TEST(RevMath, cov_sq_exp_ddv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> cov;
  std::vector<double> x(3);
  x[0] = -2;
  x[1] = -1;
  x[2] = -0.5;
  double sigma = 0.2;
      
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      stan::math::var l = 5;

      EXPECT_NO_THROW(cov = stan::math::cov_sq_exp(x, sigma, l));

      std::vector<double> grad;
      std::vector<stan::math::var> params;
      params.push_back(l);
  
      cov(i, j).grad(params, grad);

      double distance = x[i] - x[j];
      double sq_distance = stan::math::square(distance);
      double sq_l = stan::math::square(l.val());
      double exp_val = exp(sq_distance / (-2.0 * sq_l));
      EXPECT_FLOAT_EQ(stan::math::square(sigma) * exp_val, cov(i, j).val())
        << "index: (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(sigma * sigma * exp_val
                      * sq_distance/(sq_l * l.val()), grad[0]);

      stan::math::recover_memory();
    }
  }
}
