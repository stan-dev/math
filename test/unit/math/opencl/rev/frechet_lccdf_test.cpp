#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsFrechetLccdf, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, -0.6, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value(N);
  sigma_value << 0.3, 0, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value_cl(sigma_value);

  EXPECT_NO_THROW(stan::math::frechet_lccdf(y_cl, alpha_cl, sigma_cl));

  EXPECT_THROW(stan::math::frechet_lccdf(y_size_cl, alpha_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::frechet_lccdf(y_cl, alpha_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::frechet_lccdf(y_cl, alpha_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::frechet_lccdf(y_value_cl, alpha_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::frechet_lccdf(y_cl, alpha_value_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::frechet_lccdf(y_cl, alpha_cl, sigma_value_cl),
               std::domain_error);
}

auto frechet_lccdf_functor
    = [](const auto& y, const auto& alpha, const auto& sigma) {
        return stan::math::frechet_lccdf(y, alpha, sigma);
      };

TEST(ProbDistributionsFrechetLccdf, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(frechet_lccdf_functor, y, alpha,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      frechet_lccdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      sigma.transpose().eval());
}

TEST(ProbDistributionsFrechetLccdf, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd alpha(N);
  alpha << 0.5, 1.2, 1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(frechet_lccdf_functor,
                                                         y_scal, alpha, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      frechet_lccdf_functor, y_scal, alpha.transpose().eval(), sigma);
}

TEST(ProbDistributionsFrechetLccdf, opencl_broadcast_alpha) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  double alpha_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(frechet_lccdf_functor,
                                                         y, alpha_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      frechet_lccdf_functor, y.transpose().eval(), alpha_scal, sigma);
}

TEST(ProbDistributionsFrechetLccdf, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, 0.8, 1.0;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;
  double sigma_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(frechet_lccdf_functor,
                                                         y, alpha, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      frechet_lccdf_functor, y.transpose().eval(), alpha, sigma_scal);
}

TEST(ProbDistributionsFrechetLccdf, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(frechet_lccdf_functor, y, alpha,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      frechet_lccdf_functor, y.transpose().eval(), alpha.transpose().eval(),
      sigma.transpose().eval());
}
#endif
