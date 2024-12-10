#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST_F(OpenCLRevTests, probdistNormal_error_checking) {
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

  EXPECT_NO_THROW(stan::math::normal_lpdf(y_cl, mu_cl, sigma_cl));

  EXPECT_THROW(stan::math::normal_lpdf(y_size_cl, mu_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::normal_lpdf(y_cl, mu_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::normal_lpdf(y_cl, mu_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::normal_lpdf(y_value_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_lpdf(y_cl, mu_value_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::normal_lpdf(y_cl, mu_cl, sigma_value_cl),
               std::domain_error);
}

auto normal_lpdf_functor
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::normal_lpdf(y, mu, sigma);
      };
auto normal_lpdf_functor_propto
    = [](const auto& y, const auto& mu, const auto& sigma) {
        return stan::math::normal_lpdf<true>(y, mu, sigma);
      };

TEST_F(OpenCLRevTests, probdistNormal_opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(normal_lpdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(normal_lpdf_functor_propto, y,
                                                mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      normal_lpdf_functor_propto, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}

TEST_F(OpenCLRevTests, probdistNormal_opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(normal_lpdf_functor,
                                                         y_scal, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      normal_lpdf_functor_propto, y_scal, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      normal_lpdf_functor, y_scal, mu.transpose().eval(), sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      normal_lpdf_functor_propto, y_scal, mu, sigma.transpose().eval());
}

TEST_F(OpenCLRevTests, probdistNormal_opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(normal_lpdf_functor, y,
                                                         mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      normal_lpdf_functor_propto, y, mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      normal_lpdf_functor, y.transpose().eval(), mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      normal_lpdf_functor_propto, y, mu_scal, sigma.transpose().eval());
}

TEST_F(OpenCLRevTests, probdistNormal_opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  double sigma_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(normal_lpdf_functor, y,
                                                         mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      normal_lpdf_functor_propto, y, mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      normal_lpdf_functor, y.transpose().eval(), mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      normal_lpdf_functor_propto, y, mu.transpose().eval(), sigma_scal);
}

TEST_F(OpenCLRevTests, probdistNormal_opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(normal_lpdf_functor, y, mu,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(normal_lpdf_functor_propto, y,
                                                mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      normal_lpdf_functor, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      normal_lpdf_functor_propto, y.transpose().eval(), mu.transpose().eval(),
      sigma.transpose().eval());
}

#endif
