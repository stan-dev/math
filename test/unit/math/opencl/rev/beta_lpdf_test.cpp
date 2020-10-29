#ifdef STAN_OPENCL
#include <stan/math/opencl/rev/opencl.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBeta, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.1, 0.5;
  Eigen::VectorXd y_value1(N);
  y_value1 << -0.1, 0.5, 0.99;
  Eigen::VectorXd y_value2(N);
  y_value2 << 0.1, 1.1, 0.99;

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value1_cl(y_value1);
  stan::math::matrix_cl<double> y_value2_cl(y_value2);

  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);

  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::beta_lpdf(y_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::beta_lpdf(y_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::beta_lpdf(y_value1_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_value1_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_cl, beta_value1_cl),
               std::domain_error);

  EXPECT_THROW(stan::math::beta_lpdf(y_value2_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_value2_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::beta_lpdf(y_cl, alpha_cl, beta_value2_cl),
               std::domain_error);
}

auto beta_lpdf_functor
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::beta_lpdf(y, alpha, beta);
      };
auto beta_lpdf_functor_propto
    = [](const auto& y, const auto& alpha, const auto& beta) {
        return stan::math::beta_lpdf<true>(y, alpha, beta);
      };

TEST(ProbDistributionsBeta, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0, 0.5, 1;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(beta_lpdf_functor, y, alpha,
                                                beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_lpdf_functor_propto, y,
                                                alpha, beta);
}

TEST(ProbDistributionsBeta, opencl_broadcast_y) {
  int N = 3;

  double y = 0.3;
  Eigen::VectorXd y_vec(N);
  y_vec << y, y, y;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::matrix_cl<double> y_vec_cl(y_vec);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> beta_cl(beta);

  EXPECT_NEAR_REL(stan::math::beta_lpdf(y, alpha_cl, beta_cl),
                  stan::math::beta_lpdf(y_vec_cl, alpha_cl, beta_cl));
  EXPECT_NEAR_REL(stan::math::beta_lpdf<true>(y, alpha_cl, beta_cl),
                  stan::math::beta_lpdf<true>(y_vec_cl, alpha_cl, beta_cl));

  stan::math::var y_var1 = y;
  stan::math::var y_var2 = y;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var_vec(3);
  y_var_vec << y_var2, y_var2, y_var2;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha_var1 = alpha;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha_var2 = alpha;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> beta_var1 = beta;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> beta_var2 = beta;
  auto y_var_vec_cl = stan::math::to_matrix_cl(y_var_vec);
  auto alpha_var1_cl = stan::math::to_matrix_cl(alpha_var1);
  auto alpha_var2_cl = stan::math::to_matrix_cl(alpha_var2);
  auto beta_var1_cl = stan::math::to_matrix_cl(beta_var1);
  auto beta_var2_cl = stan::math::to_matrix_cl(beta_var2);

  stan::math::var res1 = stan::math::beta_lpdf(y_var1, alpha_var1_cl, beta_var1_cl);
  stan::math::var res2 = stan::math::beta_lpdf(y_var_vec_cl, alpha_var2_cl, beta_var2_cl);

  (res1 + res2).grad();

  EXPECT_NEAR_REL(res1.val(), res2.val());

  EXPECT_NEAR_REL(y_var1.adj(), y_var2.adj());
  EXPECT_NEAR_REL(alpha_var1.adj().eval(), alpha_var2.adj().eval());
  EXPECT_NEAR_REL(beta_var1.adj().eval(), beta_var2.adj().eval());
}

TEST(ProbDistributionsBeta, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double alpha = 0.3;
  Eigen::VectorXd alpha_vec(N);
  alpha_vec << alpha, alpha, alpha;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.0;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> alpha_vec_cl(alpha_vec);
  stan::math::matrix_cl<double> beta_cl(beta);

  EXPECT_NEAR_REL(stan::math::beta_lpdf(y_cl, alpha, beta_cl),
                  stan::math::beta_lpdf(y_cl, alpha_vec_cl, beta_cl));
  EXPECT_NEAR_REL(stan::math::beta_lpdf<true>(y_cl, alpha, beta_cl),
                  stan::math::beta_lpdf<true>(y_cl, alpha_vec_cl, beta_cl));

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var1 = y;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var2 = y;
  stan::math::var alpha_var1 = alpha;
  stan::math::var alpha_var2 = alpha;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha_var_vec(3);
  alpha_var_vec << alpha_var2, alpha_var2, alpha_var2;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> beta_var1 = beta;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> beta_var2 = beta;
  auto y_var1_cl = stan::math::to_matrix_cl(y_var1);
  auto y_var2_cl = stan::math::to_matrix_cl(y_var2);
  auto alpha_var_vec_cl = stan::math::to_matrix_cl(alpha_var_vec);
  auto beta_var1_cl = stan::math::to_matrix_cl(beta_var1);
  auto beta_var2_cl = stan::math::to_matrix_cl(beta_var2);

  stan::math::var res1 = stan::math::beta_lpdf(y_var1_cl, alpha_var1, beta_var1_cl);
  stan::math::var res2 = stan::math::beta_lpdf(y_var2_cl, alpha_var_vec_cl, beta_var2_cl);

  (res1 + res2).grad();

  EXPECT_NEAR_REL(res1.val(), res2.val());

  EXPECT_NEAR_REL(y_var1.adj().eval(), y_var2.adj().eval());
  EXPECT_NEAR_REL(alpha_var1.adj(), alpha_var2.adj());
  EXPECT_NEAR_REL(beta_var1.adj().eval(), beta_var2.adj().eval());
}

TEST(ProbDistributionsBeta, opencl_broadcast_beta) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  double beta = 0.3;
  Eigen::VectorXd beta_vec(N);
  beta_vec << beta, beta, beta;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> beta_vec_cl(beta_vec);

  EXPECT_NEAR_REL(stan::math::beta_lpdf(y_cl, alpha_cl, beta),
                  stan::math::beta_lpdf(y_cl, alpha_cl, beta_vec_cl));
  EXPECT_NEAR_REL(stan::math::beta_lpdf<true>(y_cl, alpha_cl, beta),
                  stan::math::beta_lpdf<true>(y_cl, alpha_cl, beta_vec_cl));

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var1 = y;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y_var2 = y;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha_var1 = alpha;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> alpha_var2 = alpha;
  stan::math::var beta_var1 = beta;
  stan::math::var beta_var2 = beta;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> beta_var_vec(3);
  beta_var_vec << beta_var2, beta_var2, beta_var2;
  auto y_var1_cl = stan::math::to_matrix_cl(y_var1);
  auto y_var2_cl = stan::math::to_matrix_cl(y_var2);
  auto alpha_var1_cl = stan::math::to_matrix_cl(alpha_var1);
  auto alpha_var2_cl = stan::math::to_matrix_cl(alpha_var2);
  auto beta_var_vec_cl = stan::math::to_matrix_cl(beta_var_vec);

  stan::math::var res1 = stan::math::beta_lpdf(y_var1_cl, alpha_var1_cl, beta_var1);
  stan::math::var res2 = stan::math::beta_lpdf(y_var2_cl, alpha_var2_cl, beta_var_vec_cl);

  (res1 + res2).grad();

  EXPECT_NEAR_REL(res1.val(), res2.val());

  EXPECT_NEAR_REL(y_var1.adj().eval(), y_var2.adj().eval());
  EXPECT_NEAR_REL(alpha_var1.adj().eval(), alpha_var2.adj().eval());
  EXPECT_NEAR_REL(beta_var1.adj(), beta_var2.adj());
}

TEST(ProbDistributionsBeta, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(beta_lpdf_functor, y, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_lpdf_functor_propto, y, alpha, beta);
}

#endif
