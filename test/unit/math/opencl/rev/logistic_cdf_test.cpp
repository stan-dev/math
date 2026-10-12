#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsLogisticCdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, -INFINITY, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value(N);
  sigma_value << 0.3, 0, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);

  EXPECT_NO_THROW(stan::math::logistic_cdf(y_cl, mu_cl, sigma_cl));

  EXPECT_THROW(stan::math::logistic_cdf(y_size_cl, mu_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::logistic_cdf(y_cl, mu_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::logistic_cdf(y_cl, mu_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::logistic_cdf(y_value_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::logistic_cdf(y_cl, mu_value_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::logistic_cdf(y_cl, mu_cl, sigma_value_cl),
               std::domain_error);
}

auto logistic_cdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::logistic_cdf(y, mu, sigma);
      };

TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      logistic_cdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}

TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_small_y_neg_inf) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, -INFINITY, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      logistic_cdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}

TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_lower_tail) {
  Eigen::VectorXd y(3);
  y << -700, -350, -1400;
  Eigen::VectorXd mu(3);
  mu << 0, 0, 0;
  Eigen::VectorXd sigma(3);
  sigma << 1, 0.5, 2;

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
}

// y == INFTY contributes a factor of 1 and zero partials; the scale partial
// used to be 0 * INFTY = NaN on the device while the CPU skipped the element.
TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_underflow_small_sigma) {
  Eigen::VectorXd y(2);
  y << 8e-298, 1.0;
  Eigen::VectorXd mu(2);
  mu << 0, 0;
  Eigen::VectorXd sigma(2);
  sigma << 1e-300, 1;

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
}

TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_y_pos_inf) {
  Eigen::VectorXd y(3);
  y << 0.3, INFINITY, 1.0;
  Eigen::VectorXd mu(3);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(3);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
}

// A scaled difference that overflows to INFTY with finite arguments
// (y = 1e10, sigma = 1e-300) contributes a factor 1 and zero partials; one
// that overflows to NEGATIVE_INFTY makes the cdf and every partial 0.
// Absolute reference values from mpmath at z = 0.5, sigma = 2, because a
// comparison with the CPU passes when both are wrong in the same way.
TEST(ProbDistributionsLogisticCdf, opencl_overflow_reference) {
  Eigen::VectorXd mu(2);
  mu << 0.0, 0.5;
  Eigen::VectorXd sigma(2);
  sigma << 1e-300, 2.0;

  for (double y0 : {1e10, -1e10}) {
    Eigen::VectorXd y(2);
    y << y0, 1.5;

    stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                  sigma);

    stan::math::var_value<stan::math::matrix_cl<double>> y_cl
        = stan::math::to_matrix_cl(y);
    stan::math::var_value<stan::math::matrix_cl<double>> mu_cl
        = stan::math::to_matrix_cl(mu);
    stan::math::var_value<stan::math::matrix_cl<double>> sigma_cl
        = stan::math::to_matrix_cl(sigma);
    stan::math::var cdf = stan::math::logistic_cdf(y_cl, mu_cl, sigma_cl);
    cdf.grad();
    Eigen::VectorXd y_adj
        = stan::math::from_matrix_cl<Eigen::VectorXd>(y_cl.adj());
    Eigen::VectorXd mu_adj
        = stan::math::from_matrix_cl<Eigen::VectorXd>(mu_cl.adj());
    Eigen::VectorXd sigma_adj
        = stan::math::from_matrix_cl<Eigen::VectorXd>(sigma_cl.adj());

    EXPECT_EQ(0.0, y_adj(0));
    EXPECT_EQ(0.0, mu_adj(0));
    EXPECT_EQ(0.0, sigma_adj(0));
    if (y0 > 0) {
      EXPECT_NEAR(0.62245933120185456, cdf.val(), 1e-15);
      EXPECT_NEAR(0.11750185610079724, y_adj(1), 1e-15);
      EXPECT_NEAR(-0.11750185610079724, mu_adj(1), 1e-15);
      EXPECT_NEAR(-0.058750928050398622, sigma_adj(1), 1e-15);
    } else {
      EXPECT_EQ(0.0, cdf.val());
      EXPECT_EQ(0.0, y_adj(1));
      EXPECT_EQ(0.0, mu_adj(1));
      EXPECT_EQ(0.0, sigma_adj(1));
    }
    stan::math::recover_memory();
  }
}

TEST(ProbDistributionsLogisticCdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(logistic_cdf_functor,
                                                         y_scal, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      logistic_cdf_functor, y_scal, mu.transpose().eval(), sigma);
}

TEST(ProbDistributionsLogisticCdf, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(logistic_cdf_functor,
                                                         y, mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      logistic_cdf_functor, y.transpose().eval(), mu_scal, sigma);
}

TEST(ProbDistributionsLogisticCdf, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double sigma_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(logistic_cdf_functor,
                                                         y, mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      logistic_cdf_functor, y.transpose().eval(), mu, sigma_scal);
}

TEST(ProbDistributionsLogisticCdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(logistic_cdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      logistic_cdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}
#endif
